#!/usr/bin/env python
# Archiving script Version 3.1

# stdlib
import os
import pathlib
import subprocess
import termios
import sys
import tty
import re
import datetime
import time
import hashlib
import sqlite3
import getpass
import pickle
import pwd
import grp
import shutil
import logging
from itertools import compress, cycle

# packages
import pytz
import click
import pandas as pd


logtime = datetime.datetime.now().strftime('%d-%b-%Y_%H.%M')

__version__ = '4.4'

def get_type(ext, size, settings):
    types = settings['types']
    sizes = settings['sizes']

    # Turn types inside out for access
    extensions = {}
    for k, v in types.items():
        for vv in v:
            extensions[vv] = k

    # Determine type
    if ext.lower() in extensions.keys():
        kind = extensions[ext.lower()]
        if size > sizes[kind]:
            kind = 'other'
    else:
        kind = 'other'

    return kind


def getinfo_dir(entry):
    keep = 'unknown'
    if type(entry) == os.DirEntry:
        path = entry.path
    else:
        path = str(entry)
    ext = 'directory'
    directory = True
    if entry.name == '.git':
        stats = entry.stat()
        mode = stats.st_mode
        size_mib = stats.st_size / 1048576 if not entry.is_symlink() else 0
        keep = 'git'
        kind = 'gitdir'
        modified = None
        created = None
        user = None
        group = None
    elif entry.is_symlink():
        return None
    else:
        keep = 'directory'
        kind = 'directory'
        size_mib = -9
        mode = -9
        modified = None
        created = None
        user = None
        group = None
    return {'size_mib': size_mib, 'mode': mode, 'path': path, 'keep': keep,
            'filename': entry.name, 'extension': ext, 'kind': kind,
            'directory': directory, 'time_modified': modified, 'user': user,
            'group': group}


def getinfo_file(entry, settings):
    if type(entry) == os.DirEntry:
        path = entry.path
    else:
        path = str(entry)
    keep = 'unknown'

    if entry.name == 'Snakefile':
        ext = '.smk'
    elif entry.name == '.gitignore':
        ext = '.gitignore'
    else:
        ext = os.path.splitext(entry.name)[1]

    directory = False
    pathparts = pathlib.Path(path).parts
    pardir = '.' if len(pathparts) == 1 else pathparts[-2]

    if entry.is_symlink() and not os.path.exists(path):
        kind = 'broken_link'
        keep = 'no'
        size_mib = 0
        mode = None
        modified = None
        created = None
        user = None
        group = None
    elif not os.path.exists(path):
        kind = 'fs_err'
        keep = 'yes'
        size_mib = 0
        mode = None
        modified = None
        created = None
        user = None
        group = None
    else:
        stats = entry.stat()
        mode = stats.st_mode
        size_mib = stats.st_size / 1048576 if not entry.is_symlink() else 0
        modified = stats.st_mtime
        created = stats.st_ctime
        user = stats.st_uid
        group = stats.st_gid
        if '.snakemake' in pathparts:
            if pardir == 'log':
                keep = 'snakemake_log' if pardir == 'log' else 'no'
                kind = 'log'
            else:
                keep = 'no'
                kind = 'other'
        elif pardir == '.snakejob':
            keep = 'snakejob_log'
            kind = 'log'
        elif '.git' in pathparts or ext == '.gitignore':
            keep = 'git_file'
            kind = 'gitfile'
        else:
            kind = get_type(ext, size_mib, settings)
            if kind == 'other' and size_mib > 0:
                keep = 'no'
            elif size_mib == 0:
                keep = 'empty'
            else:
                keep = kind
    return {'size_mib': size_mib, 'mode': mode, 'path': path, 'keep': keep,
            'filename': entry.name, 'extension': ext, 'kind': kind,
            'directory': directory, 'time_modified': modified, 'user': user,
            'group': group}


def getinfo(entry, settings):
    try:
        if entry.is_dir():
            return getinfo_dir(entry)
        else:
            return getinfo_file(entry, settings)
    except PermissionError:
        if not entry.is_symlink():
            raise
        return {'size_mib': 0, 'mode': None, 'path': entry.path, 'keep': 'no',
                'filename': entry.name, 'extension': None, 'group': None,
                'kind': 'forbidden_link', 'directory': None, 'user': None,
                'time_modified': None}


def git_test(gitdir):
    with open(gitdir + '/info/exclude', 'a') as excludefile:
        print('archive.log', file=excludefile)
    git_status = subprocess.run(['git -c color.status=always status',
                                 gitdir + '/..'],
                                capture_output=True, shell=True)
    gs_return = git_status.stdout.decode()
    utd = 'Your branch is up to date' in gs_return
    utd_or_behind = utd or 'Your branch is behind' in gs_return
    ahead = 'Your branch is ahead' in gs_return
    changes = 'Changes not staged' in gs_return
    untracked = 'Untracked files' in gs_return

    if ahead or changes or untracked:
        git_issues = list(compress(
            ['unpushed commits', 'unstaged chages', 'untracked files'],
            [ahead, changes, untracked]))
        if len(git_issues) == 1:
            git_issues = git_issues[0]
        else:
            git_issues = '{} and {}'.format(
                ', '.join(git_issues[:-1]),  git_issues[-1])
        msg = ('\033[93mThe git repo \033[1m{repo}\033[22m'
               + ' has {iss}\033[0m:\n\n')
        print(msg.format(repo=gitdir, iss=git_issues)
              + gs_return)
        return 2 if changes or untracked else 1
    return 0


def git_test_and_prompt(gitdir, batch=False):
    git_status = git_test(gitdir)
    if git_status == 1:
        print("We recommend you cancel the script and commit changes.")
        if not batch and question("Do you want to cancel?", True):
            logging.info(
                f'{gitdir} is a git project and is ahead of remote. Exiting')
            exit(1)
        else:
            logging.warning(
                f'{gitdir} is a git project and is ahead of remote')
    elif git_status == 2:
        if not batch and question("Do you want to cancel to commit and push changes?", True):
            logging.warning(
                f'{gitdir} is a git project and has untracked changes')
            exit(1)
        else:
            logging.warning(
                f'{gitdir} is a git project and has untracked changes')
    else:
        logging.info(f'{gitdir} is a git project and is up to date')


def question_old(text, default=None, prepend=''):
    '''prompt with default. Use click.confirm if converting to package'''
    def _getch(message):
        print(message)
        fd = sys.stdin.fileno()
        old_settings = termios.tcgetattr(fd)
        try:
            tty.setraw(fd)
            ch = sys.stdin.read(1)     # This number represents the length
        finally:
            termios.tcsetattr(fd, termios.TCSADRAIN, old_settings)
        return ch

    if isinstance(default, bool):
        yn = ' [Y/n]' if default else ' [y/N]'
    else:
        yn = ' [y/n]'

    try:
        answer = _getch(prepend + text + yn)
    except:
        answer = input(prepend + text + yn)
    if answer in ['', '\r'] and isinstance(default, bool):
        return default
    elif answer in ['', '\r']:
        question_old(text, default, 'Please enter a response:\n')
    elif answer.lower()[0] in ["y", "yes"]:
        return True
    elif answer.lower()[0] in ["n", "no"]:
        return False
    else:
        question_old(text, default, 'Invalid answer; please try again:\n')

def question(text, default=None, prepend='', batch=False, default_batch=None):
    if batch:
        if default_batch is None and default in [True, False]:
            return default
        elif default_batch in [True, False]:
            return default_batch
        else:
            logging.error(f"In batch mode with no default for question {text}")
            raise ValueError(f"In batch mode with no default for question {text}")
    return click.confirm(prepend + text, default=default)


def load_config(configfile):
    if os.path.isfile(configfile):
        import yaml
        with open(configfile, 'r') as conf:
            settings = yaml.safe_load(conf)
    else:
        settings = {
            'types': {
                'log': ['.log', '.err', '.stderr', '.out', '.stdout'],
                'readme': ['.readme'],
                'config': ['.yaml', '.yml', '.conf', '.json', '.ini'],
                'script': ['.sh', '.r', '.rmd', '.smk', '.py'],
                'txt_and_out': ['.txt', '.pdf', '.html',
                                '.htm', '.org', '.md'],
                'image': ['.tiff', '.png', '.jpg', '.svg']
                },
            'sizes': {
                'log': 200, 'readme': 10, 'config': 50, 'script': 5,
                'txt_and_out': 100, 'image': 5
                }
            }
    return settings


def scantree(path, settings=load_config(False)):
    '''Recursively yield DirEntry objects for given directory.'''
    for entry in os.scandir(path):
        if entry.is_dir(follow_symlinks=False):
            if entry.name == '.git':
                yield getinfo(entry, settings)
            yield from scantree(entry.path, settings)
        else:
            yield getinfo(entry, settings)


def scanfiles(paths='.', settings=load_config(False)):
    '''Recursively yield DirEntry objects for given directories and files.'''
    if type(paths) is str:
        paths = [paths]
    for path in paths:
        path = pathlib.Path(path)
        if path.is_dir():
            if path == '.git':
                yield getinfo(path, settings)
            yield from scantree(path, settings)
        else:
            yield getinfo(path, settings)

def list_files(path='.', settings=load_config(False)):
    files = [x for x in scanfiles(path, settings) if x is not None]
    files = pd.DataFrame.from_records(files)
    files['user'] = [pwd.getpwuid(int(x))[0]
                     if not pd.isnull(x) else None
                     for x in files['user']]
    files['group'] = [grp.getgrgid(int(x))[0]
                      if not pd.isnull(x) else None
                      for x in files['group']]
    files['time_modified'] = [datetime.datetime.fromtimestamp(x)
                              if not pd.isnull(x) else None
                              for x in files['time_modified']]
    return files


def check_tmux():
    '''Check that screen or tmux are in use:'''
    if os.getenv('TERM') != 'screen':
        print('We recommend running in screen or tmux.')
        logging.warning('Screen/TMUX session not detected.')
        if question('Quit to run screen or tmux?', True):
            logging.error('Quitting because not in Screen or TMUX.')
            exit(2)


def get_sizes(files):
    files_nobroken = files[files['mode'].notnull()]
    files_nobroken = files_nobroken.astype({'mode': 'int32'})

    total = files.loc[:, 'size_mib'].sum().round(1)
    freed = files[files['keep'] == 'no'].loc[:, 'size_mib'].sum().round(1)
    kept = (total - freed).round(1)

    kinds = files_nobroken.groupby(['kind']).sum(numeric_only=True).round(3)[['size_mib']]

    print('The total space taken up by each kind of file is as follows:\n')
    print(kinds)
    print('\n\nThe total size of your directory is {} MiB.'.format(int(total)))
    print('''If you choose to keep small files, you will free at least {} MiB
      leaving at most {} MiB.'''.format(int(freed), int(kept)))

    logging_usage = 'File usage by type:\n\n{}'.format(kinds)

    sizes = {'total': total, 'freed': freed, 'kept': kept}
    return logging_usage, sizes


def test_tar():
    '''
    Test to see if the available tar is new enough to suport multithreading
    and ACLs.
    '''
    output = subprocess.run(['tar', '--version'],
                            capture_output=True).stdout
    vs = output.splitlines()[0].decode()
    pat = r'(\d+)\.([0-9.]+)$'
    ver = {}
    ver['maj'], ver['min'] = re.search(pat, vs).groups()
    ver = {x : float(y) for x, y in ver.items()}
    if ver['maj'] > 1 or ver['min'] > 28:
        return True
    return False


def exit_code(logname):
    with open(logname, 'r') as f:
        for line in f.readlines():
            ecodematch = re.search(r'Exited with exit code (\d+)\.', line)
            if ecodematch:
                return int(ecodematch.groups()[0])
    return 0


def make_tarball(files, total, batch=False, wdir='.',
                 date_time=None, fname_xtra=None,
                 indiv_files=False, compression=True):
    if wdir == '.':
        wdir = os.getcwd()
    fullpath = os.path.normpath(wdir)
    parent = os.path.dirname(fullpath)
    cwd = os.path.basename(fullpath)
    try:
        tznm = '/'.join(os.readlink('/etc/localtime').split('/')[-2:])
    except:
        tznm = 'America/New_York'
    dt_current = datetime.datetime.now(tz=pytz.timezone(tznm))
    date = dt_current.strftime('%G-%b-%d')
    dt_iso = dt_current.isoformat()
    hash_dtd = filehash(cwd, dt_iso)
    files['file_id'] = [filehash(x, dt_iso) for x in files['path']]
    files['archive_id'] = hash_dtd
    dt = date_time if date_time else date
    dte = f'{fname_xtra}_{dt}' if fname_xtra else dt
    tarball = f'{cwd}_{dte}.tar.bz2' if compression else f'{cwd}_{dte}.tar'
    username = getpass.getuser()

    archive_info = {
        'archive_id': hash_dtd,
        'user_name': username,
        'time': dt_current,
        'archive_name': tarball,
        'archive_directory': fullpath,
        'total_mib': total
        }


    jobname = f'{cwd}_{dte}.tarball'
    oo = f'{jobname}.stdout'
    eo = f'{jobname}.stderr'

    if indiv_files:
        file_list = [os.path.join(cwd, x) for x in files['path']]
        tar_list = ' '.join(file_list)
    else:
        tar_list = cwd

    tempdir = '/sc/arion/scratch/{}/archiving'.format(username)

    newtar = test_tar()

    temp_tarball = '{}/{}'.format(tempdir, tarball)
    if compression:
        tar_cmd = 'tar -cvpjf {} -C {} {}'.format(temp_tarball, parent, tar_list)
        tar_cmd_fast = 'tar c -I"pbzip2 -p24" -f {} -C {} --xattrs {}'
        tar_cmd_fast = tar_cmd_fast.format(temp_tarball, parent, tar_list)
    elif newtar:
        tar_cmd = f'tar -cvpf {temp_tarball} -C {parent} --xattrs {tar_list}'
    else:
        tar_cmd = f'tar -cvpf {temp_tarball} -C {parent} {tar_list}'

    try:
        not_used = subprocess.run(["pbzip2", "-h"], capture_output=True)
    except FileNotFoundError:
        pbzip = False
    else:
        pbzip = True


    ncore = 1
    if not compression:
        logging.info("Not compressing tarball.")
        if newtar:
            logging.info("Using tar to archive with XATTRs.")
    elif newtar and pbzip:
        logging.info("Using pbzip2 to compress in parallel. XATTRs recorded.")
        tar_cmd = tar_cmd_fast
        ncore = 24
    elif pbzip:
        print("Tar is older than v1.28, so multithreading tar_list unsupported.")
        print("Use 'mamba install -c conda-forge tar' to enable features.")
        logging.warning("Tar is older than 1.28. Using serial compression")
    elif newtar:
        print("Parallel bzip2 missing. Archiving will be slow unless you install")
        print("pbzip2. Use 'mamba install -c conda-forge pbzip2' to enable.")
        logging.warning("pbzip2 is missing. Using slow serial compression.")
    else:
        print("Parallel bzip2 missing and tar is older than v1.28 so multithreading")
        print("is unsupported.")
        print("Use 'mamba install -c conda-forge pbzip2 tar' to enable features.")
        logging.warning("pbzip2 is missing and tar is older than 1.28.")
        logging.warning("Using slow serial compression.")

    job_cmd = ('bsub -P acc_LOAD -q premium -n {} -R rusage[mem=1000] '
               + '''-R span[hosts=1] -W 140:00 -J {} -oo {} -eo {} '{}' ''')
    job_cmd = job_cmd.format(ncore, jobname, oo, eo, tar_cmd)
    tempdir_path = pathlib.Path(tempdir)
    pathlib.Path.mkdir(tempdir_path, parents=True, exist_ok=True)
    tar_incomplete = True
    tar_attempt = 1

    while tar_incomplete and tar_attempt < 4:
        archive_job = subprocess.run(job_cmd, capture_output=True, shell=True)
        print(archive_job.stdout.decode())
        print(archive_job.stderr.decode())
        archive_job.check_returncode()
        job_expression = r'(?<=^Job <)\d+(?=> is submitted to queue)'
        job = re.search(job_expression, archive_job.stdout.decode())[0]
        completed = False
        spinner = cycle(['◐', '◓', '◑', '◒'])
        print("Waiting to check job.", end='\r')
        while not completed:
            time.sleep(5)
            jobstat = subprocess.run(
                'bjobs -do "stat" -noheader {}'.format(job),
                capture_output=True, shell=True)
            jobstat = jobstat.stdout.decode().strip()
            if jobstat == 'PEND':
                statmsg = "Tar job is pending... {}".format(next(spinner))
                print(statmsg, end='\r')
            elif jobstat == 'RUN':
                statmsg = "Tar job is running... {}".format(next(spinner))
                print(statmsg, end='\r')
            else:
                statmsg = "Tar job has completed.   "
                print(statmsg, end='\n')
            completed = jobstat not in ['RUN', 'PEND']

        tar_code = exit_code(oo)
        if jobstat != 'EXIT' and tar_code == 0:
            tar_incomplete = False
            break
        elif jobstat == 'EXIT':
            logging.warning(
                'LSF archive job failed on attempt {}'.format(tar_attempt))
            print('LSF archive job failed on attempt {}'.format(tar_attempt))
            if tar_attempt < 3:
                print('Trying again up to {} times'.format(3 - tar_attempt))
            else:
                logging.error('Aborting because LSF tar job failed.')
                raise Exception('Aborting because LSF tar job failed.')
            tar_attempt += 1
        else:
            logging.warning(
                'tar exited with exit code {} on attempt {}'.format(
                    tar_code, tar_attempt))
            print('tar exited with exit code {} on attempt {}'.format(
                    tar_code, tar_attempt))
            if tar_code == 1:
                print('The following files were unsuccessful:')
                with open(eo, 'r') as f:
                    print(f.read(), end='\n')
                if not batch:
                    print('You can choose to try again, to use this archive,'
                          + '\n or to abort all operations.')
                if tar_attempt < 3 and question('Try again?', True, batch=batch):
                    print('Trying archiving again.')
                    logging.warning('Trying archiving again.')
                    tar_attempt += 1
                elif batch:
                    logging.error('Aborting because tar is incomplete in batch mode.')
                    os.remove(temp_tarball)
                    raise Exception('Aborting because tar is incomplete.')
                elif question('Continue with incomplete archive?'):
                    logging.warning('Continuing with incomplete archive.')
                    print('WARNING: Continuing with incomplete archive.')
                    break
                else:
                    logging.error('Aborting because tar is incomplete.')
                    os.remove(temp_tarball)
                    raise Exception('Aborting because tar is incomplete.')
            else:
                retry = tar_attempt < 3
                if batch:
                    logging.warning('Archiving failed.')
                else:
                    print('Archiving failed. You can choose to try again'
                          + '\n or to abort all operations.')
                    retry = retry and question('Try again?', True)
                if retry:
                    print('Trying archiving again.')
                    logging.warning('Trying archiving again.')
                else:
                    logging.error('Aborting because tar failed.')
                    raise Exception('Aborting because tar failed.')
            tar_attempt += 1
    return archive_info, temp_tarball


def filehash(file, date):
    '''Generate hash from file/dir name and date'''
    hash = hashlib.shake_256()
    hash.update((file + date).encode())
    return hash.hexdigest(8)


def file_rm(files, archive_info, sizes, log_sz, keep, logname=None):
    logging.info('Removing files')
    files['removal'] = 'not removed'
    if keep in ['ask', 'default']:
        keep_tf = question('Do you want to keep {} MiB of small files?'.format(
                           sizes['kept']), True)
    elif keep == 'none':
        keep_tf = False
    else:
        keep_tf = True
    if keep_tf:
        logging.info('Keeping {} MiB of small files'.format(sizes['kept']))
        logprint(log_sz, logname)
        archive_info['freed_mib'] = sizes['freed']
        archive_info['kept_mib'] = sizes['kept']
        for file in files[files['keep'] == 'no']['path']:
            try:
                os.remove(file)
                files.loc[files['path'] == file, 'removal'] = 'success'
            except IsADirectoryError:
                files.loc[files['path'] == file, 'removal'] = 'fail: dir'
            except PermissionError:
                files.loc[files['path'] == file, 'removal'] = 'fail: perm'
            except FileNotFoundError:
                files.loc[files['path'] == file, 'removal'] = 'fail: miss'
            except OSError:
                files.loc[files['path'] == file, 'removal'] = 'fail: OS'
        files.loc[files['keep'] != 'no', 'removal'] = 'kept'
    else:
        logging.info('Not keeping small files'.format(sizes['kept']))
        archive_info['freed_mib'] = sizes['total']
        archive_info['kept_mib'] = 0
        for file in files['path']:
            try:
                os.remove(file)
                files.loc[files['path'] == file, 'removal'] = 'success'
            except IsADirectoryError:
                files.loc[files['path'] == file, 'removal'] = 'fail: dir'
            except PermissionError:
                files.loc[files['path'] == file, 'removal'] = 'fail: perm'
            except FileNotFoundError:
                files.loc[files['path'] == file, 'removal'] = 'fail: miss'
            except OSError:
                files.loc[files['path'] == file, 'removal'] = 'fail: OS'
    logging.info('Done removing files')
    return files, archive_info


def delete_files(files, archive_info, sizes, log_sz, keep,
                 batch=False, logname=None):
    idx = files[files['filename'] == 'archive_{}.log'.format(logtime)].index
    files.drop(idx, inplace=True)
    if batch or question('Do you still want to remove files after generating archive'):
        try:
            files, archive_info = file_rm(files, archive_info, sizes, log_sz,
                                          keep, logname)
        except:
            logging.exception('Uncaught exception deleting files')
            print("Unexpected error deleting files:", sys.exc_info()[0])
            with open('archive_dump.p', 'wb') as pklh:
                pickle.dump(files, pklh)
            if batch or question('Do you want to continue anyway?', False):
                files['removal'] = 'fail: other'
                archive_info['freed_mib'] = 0
                archive_info['kept_mib'] = sizes['total']
                archive_info['removal_failure'] = 'yes'
            else:
                archive_info['removal_failure'] = 'yes'
                with open('info_dump.p', 'wb') as pklh:
                    pickle.dump(archive_info, pklh)
                raise
        removed = files.loc[files['removal'] == 'success', 'size_mib'].sum()
        removed = removed.round()
    else:
        logging.info('Not removing files')
        archive_info['freed_mib'] = 0
        archive_info['kept_mib'] = sizes['total']
        files['removal'] = 'kept'
        removed = 0
    discrep = round(archive_info['freed_mib']) - removed
    removed = 'Removed {} MiB of files.'.format(int(removed))
    print(removed)
    logging.info(removed)
    if discrep != 0:
        discrep = 'Could not remove {} MiB of files.'.format(int(discrep))
        print(discrep)
        logging.warning(discrep)
    return files, archive_info

def delete_no_files(files, archive_info, sizes):
    idx = files[files['filename'] == 'archive_{}.log'.format(logtime)].index
    files.drop(idx, inplace=True)
    logging.info('Not removing files')
    archive_info['freed_mib'] = 0
    archive_info['kept_mib'] = sizes['total']
    files['removal'] = 'kept'
    return files, archive_info


def tsm_archive(temp_tarball, archive_info, files=None, attempt=1,
                keep_tar=False, batch=False, ans_codes=[], logname=None):
    archive_name = archive_info['archive_name']
    if attempt == 1:
        logging.info('Moving archive tarball')
        print('Moving archive tarball')
        shutil.copy(temp_tarball, archive_name)
        os.remove(temp_tarball)
    logging.info('DSMC Job starting')
    print('DSMC Job starting')
    dsmc_cmd = [shutil.which('dsmc'), 'archive',
                '-se={}'.format(archive_info['user_name']),
                '"{}"'.format(archive_name)]
    with subprocess.Popen(dsmc_cmd, stdout=subprocess.PIPE,
                          stderr=subprocess.STDOUT,
                          bufsize=1, universal_newlines=True,
                          ) as dsmc_job:
        for line in dsmc_job.stdout:
            print(line, end='')
            logprint(line.strip(), logname)
            ans_codes += re.findall(r'ANS\d{4}[WES]', line)
    logging.info('DSMC Job completed')
    warnings = {
        'OK': {
            'ANS1144W': 'No password auth',
            'ANS1261W': 'archive description is empty string',
            'ANS1424W': 'Retrying after issue',
            'ANS1587W': 'Unable to read xattrs for file',
            'ANS1740W': 'Unable to read ACLs for file',
            'ANS1741W': 'Unable to read xattrs for file',
            'ANS1809W': 'TSM disconnect'
            },
        'retry': {
            'ANS1588W': 'IO error reading attrs'
            },
        'fail': {
            'ANS1133W': "Bad wildcard",
            'ANS1184W': "Command not supported",
            'ANS1115W': 'Filename excluded',
            'ANS1252W': 'Function not supported on server',
            'ANS1375W': 'User requested skip', # We did not
            'ANS1832W': 'Option no longer supported',
            'ANS1946W': 'File exists; skipping',
            'ANS1947W': 'Directory exists; skipping'
        }
    }
    if dsmc_job.returncode == 8:
        logging.warning(f'DSMC job (attempt {attempt}) completed with warnings.')
        if batch and attempt < 4:
            print("DSMC Job (attempt {attempt}) had warnings with return {dsmc_job.returncode}. Trying again.")
            logging.info("Trying again.")
            tsm_archive(temp_tarball, archive_info, files,
                        attempt + 1, keep_tar, batch=True,
                        ans_codes=ans_codes, logname=logname)
        elif attempt < 4 and question(f'Archiving failed with return {dsmc_job.returncode}. Try again?', False):
            tsm_archive(temp_tarball, archive_info, files,
                        attempt + 1, keep_tar, ans_codes=ans_codes,
                        logname=logname)
        else:
            ans_codes = list(set(ans_codes))
            ans_codes_n = len(ans_codes)
            ans_codes_s = 's' if ans_codes_n > 1 else ''
            ans_code_waswere = 'was' if ans_codes_n == 1 else 'were'
            ans_codes = ', '.join(ans_codes[:-1]) + ' and ' + ans_codes[-1] if ans_codes_n > 1 else ans_codes[0]
            ans_code_msg = f"The TSM ANS code{ans_codes_s} {ans_code_waswere} {ans_codes}."
            if files is not None:
                with open('archive_dump.p', 'wb') as pklh:
                    pickle.dump(files, pklh)
            archive_info['exception'] = 'archiving'
            with open('info_dump.p', 'wb') as pklh:
                pickle.dump(archive_info, pklh)
            logging.exception('Uncaught warning archiving files')
            logging.exception(ans_code_msg)
            logging.exception('Please contact Brian to continue the archiving.')
            print("Unexpected warning(s) archiving files:", sys.exc_info()[0])
            print(ans_code_msg)
            print("\033[1m\033[91mContact Brian to finish archiving!\033[0m")
            raise subprocess.CalledProcessError(dsmc_job.returncode, dsmc_cmd)
    elif dsmc_job.returncode != 0:
        logging.warning(f"DSMC Job (attempt {attempt}) was unsuccessful with return {dsmc_job.returncode}.")
        if batch and attempt < 4:
            print("DSMC Job (attempt {attempt}) was unsuccessful with return {dsmc_job.returncode}. Trying again.")
            logging.info("Trying again.")
            tsm_archive(temp_tarball, archive_info, files,
                        attempt + 1, keep_tar, batch=True,
                        ans_codes=ans_codes, logname=logname)
        elif attempt < 4 and question(f'Archiving failed with return {dsmc_job.returncode}. Try again?', False):
            tsm_archive(temp_tarball, archive_info, files,
                        attempt + 1, keep_tar, ans_codes=ans_codes,
                        logname=logname)
        else:
            ans_codes = list(set(ans_codes))
            ans_codes_n = len(ans_codes)
            ans_codes_s = 's' if ans_codes_n > 1 else ''
            ans_code_waswere = 'was' if ans_codes_n == 1 else 'were'
            # join the set of ans_codes with ", " for for all but the last, and " and " for the last
            ans_codes = ', '.join(ans_codes[:-1]) + ' and ' + ans_codes[-1] if ans_codes_n > 1 else ans_codes[0]
            ans_code_msg = f"The TSM ANS code{ans_codes_s} {ans_code_waswere} {ans_codes}."
            if files is not None:
                with open('archive_dump.p', 'wb') as pklh:
                    pickle.dump(files, pklh)
            archive_info['exception'] = 'archiving'
            with open('info_dump.p', 'wb') as pklh:
                pickle.dump(archive_info, pklh)
            logging.exception('Uncaught exception archiving files')
            logging.exception(ans_code_msg)
            logging.exception('Please contact Brian to continue the archiving.')
            print("Unexpected error archiving files:", sys.exc_info()[0])
            print(ans_code_msg)
            print("\033[1m\033[91mContact Brian to finish archiving!\033[0m")
            raise subprocess.CalledProcessError(dsmc_job.returncode, dsmc_cmd)
    elif keep_tar:
        logging.info('DSMC Job successful; keeping tarball.')
    else:
        logging.info('DSMC Job successful; deleting tarball.')
        os.remove(archive_name)


def write_database(dbpath, files, archive_info):
    try:
        if 'removal_failure' in archive_info:
            del archive_info['removal_failure']
        logging.info('Writing to lab archive DB')
        engine = sqlite3.connect(dbpath)
        files.to_sql('file', con=engine, if_exists="append", index=False)
        archive_info = pd.DataFrame(archive_info, [0])
        archive_info.to_sql('archive', con=engine,
                            if_exists="append", index=False)
    except:
        with open('archive_dump.p', 'wb') as pklh:
            pickle.dump(files, pklh)
        archive_info['exception'] = 'db'
        with open('info_dump.p', 'wb') as pklh:
            pickle.dump(archive_info, pklh)
        logging.exception('Exception writing to lab DB')
        print("Unexpected error archiving files:", sys.exc_info()[0])
        print("Please contact Brian at brian.fulton-howard@mssm.edu!")
    else:
        logging.info('Done writing to lab archive DB')


def write_tables(files, archive_info, fname_extra=None, logname=None):
    try:
        if type(fname_extra) is str:
            fname_info = f'archive_info_{fname_extra}_{logtime}.log'
            fname = f'archived_file-list_{fname_extra}_{logtime}.tsv.gz'
        else:
            fname_info = f'archive_info_{logtime}.log'
            fname = f'archived_file-list_{logtime}.tsv.gz'
        logging.info('Writing text tables to directory')
        files.to_csv(fname, sep='\t', index=False)
        archive_info = '''
Archive ID: {archive_id},
Username: {user_name},
Archive creation time: {time},
Archive filename: {archive_name},
Original Directory: {archive_directory},
Total original size (MiB): {total_mib}
'''.format(**archive_info)
        print(archive_info)
        logprint(archive_info, logname)
        with open(fname_info, 'w') as ai:
            print(archive_info, file=ai)
    except:
        with open('archive_dump.p', 'wb') as pklh:
            pickle.dump(files, pklh)
        archive_info['exception'] = 'tables'
        with open('info_dump.p', 'wb') as pklh:
            pickle.dump(archive_info, pklh)
        logging.exception('Exception writing text tables to directory')
        print("Unexpected error outputing tables", sys.exc_info()[0])
        print("Please contact Brian at brian.fulton-howard@mssm.edu!")
        raise
    else:
        logging.info('Done writing text tables to directory')


def logprint(message, logname):
    logname = f'archive_{logtime}.log' if logname is None else logname
    with open(logname, 'a') as logfile:
        print('\n{}\n'.format(message), file=logfile)

main_opts = {
    'd': click.option('-d', '--delete', default=False, is_flag=True,
        help='Delete files after confirming'),
    'k': click.option('-k', '--keep', 
        type=click.Choice(['small', 'none', 'ask', 'default'],
                          case_sensitive=False),
        default='default',
        help='Keep small files, none, or ask (default is ask when deleting)'),
    'K': click.option('-K', '--keep-config', type=click.File('r'),
        help='Use YAML file to determine which files to keep'),
    't': click.option('-t', '--keep-tarball',
        type=click.Choice(['yes', 'no', 'ask'], case_sensitive=False),
        default='ask',
        help='Keep tarball after archiving'),
    'u': click.option('-u', '--uncompressed', is_flag=True, default=False,
        help='Do not compress tarball with bzip2'),
    'f': click.option('-f', '--files', multiple=True,
        help='Files to archive'),
    }

def main_opt_get(k):
    return main_opts[k]

# Template for all command line usage
def archive(delete, keep="ask", keep_config=None, keep_tarball="no",
            safe=False, batch=False, directory=".", files=None,
            compression=True):
    files = files if files else None
    if directory != ".":
        # get current directory
        dir_run = os.getcwd()
        # change to the directory
        os.chdir(directory)
        # get working directory
        dir_archive = os.getcwd()
    else:
        dir_run = os.getcwd()
        dir_archive = dir_run
    
    if batch:
        logtime = datetime.datetime.now().strftime('%d-%b-%Y_%Hh%Mm%Ss%fus')

    logname = f'archive_{logtime}.log'

    logging.basicConfig(
        filename=logname, level=logging.INFO,
        format='%(asctime)s %(levelname)s: %(message)s',
        datefmt='%m/%d/%Y %H:%M:%S')

    logging.info(f"Running in {dir_run}")
    logging.info(f"Archiving in {dir_archive}")
    print(f"Achiving in {dir_archive}")

    if not batch:
        check_tmux()

    if not delete and keep != 'default':
        print('You must specify "-d" or "--delete" to delete files')
        logging.critical(
            '"-d" or "--delete" not specified but keep is not default')
        exit(1)
    if not delete and keep_config is not None:
        print('You must specify "-d" or "--delete" to delete files')
        logging.critical(
            '"-d" or "--delete" not specified but keep file specified')
        exit(1)
    if keep_config is not None and keep != 'default':
        print('You cannot specify both "-k" and "-K"')
        logging.critical('Both "-k" and "-K" specified')
        exit(1)
    if safe and delete:
        print('You cannot delete files in safe mode.')
        logging.critical('"-d" or "--delete" specified in safe mode.')
        exit(1)
    if batch and keep_tarball == 'ask':
        logging.critical('"keep_tarball" set to "ask" in batch mode.')
        raise ValueError('You cannot ask about tarball keeping in batch mode.')
    if batch and delete and keep in ['ask', 'default']:
        print('You cannot ask about tarball keeping in batch mode.')
        logging.critical('"keep" set to "ask" or "default" in batch mode.')
        raise ValueError('You cannot ask about keeping files in batch delete mode. Valid options for "keep" are "small" or "none"')

    if type(keep_config) is dict:
        settings = keep_config
    elif keep_config:
        settings = load_config(keep_config)
    else:
        settings = load_config(False)

    print('Scanning for files')
    logging.info(f'Archiving script v{__version__}')
    logging.info('Scanning directory for files')
    if files is None:
        fname_extra = None
        files = list_files(settings)
        individual_files = False
    else:
        if type(files) is str:
            fname_extra = files
            files = [files]
        elif type(files) is list:
            if len(files) < 6:
                fname_extra = '_-_'.join(os.path.basename(x) for x in files)
            else:
                fname_extra = f'{len(files)}-files'
        else:
            logging.error('Files argument must be a string or list of strings')
            raise ValueError('Files argument must be a string or list of strings')
        for file in files:
            if not os.path.exists(file):
                logging.error(f'File {file} does not exist.')
                raise FileNotFoundError(f'File {file} does not exist.')
            if file[0] == '/':
                logging.error(f'File {file} is not a relative path.')
                raise ValueError(f'File {file} is not a relative path.')
        files = ["./" + str(pathlib.PurePath(file)) for file in files]
        files = list_files(files, settings)
        individual_files = True
    logging.info('Done scanning directory for files')

    if any(files['kind'] == 'gitdir') and delete:
        logging.info("Git directories detected")
        for gitdir in files[files['kind'] == 'gitdir']['path']:
            git_test_and_prompt(gitdir, batch)

    log_sz, sizes = get_sizes(files)

    if batch or question('Archive this folder:\n{}?'.format(os.getcwd())):
        archive_info, temp_tarball = make_tarball(files, sizes['total'], batch,
                                                  date_time=logtime,
                                                  fname_xtra=fname_extra,
                                                  wdir=dir_archive,
                                                  indiv_files=individual_files,
                                                  compression=compression)
        if delete == True:
            files, archive_info = delete_files(files, archive_info, sizes,
                                               log_sz, keep, batch=batch,
                                               logname=logname)
        else:
            files, archive_info = delete_no_files(files, archive_info, sizes)
        if keep_tarball == 'ask' and delete == True:
            keep_tarball_tf = question(
                'Keep archive tarball in directory after TSMC archiving?' +
                '(may not reduce space used by directory)',
                default=False)
        elif keep_tarball == 'ask':
            keep_tarball_tf = question(
                'Keep archive tarball in directory after TSMC archiving?' +
                '(NOT RECOMMENDED: may double disk usage)',
                default=False)
        elif keep_tarball == 'yes':
            keep_tarball_tf = True
        else:
            keep_tarball_tf = False
        tsm_archive(temp_tarball, archive_info, files, keep_tar=keep_tarball_tf,
                    batch=batch, logname=logname)
        database = '/sc/arion/projects/LOAD/archive/archive.sqlite'
        write_database(database, files, archive_info)
        write_tables(files, archive_info, fname_extra, logname)
        logging.info('Done')
        if 'removal_failure' in archive_info:
            with open('info_dump.p', 'wb') as pklh:
                pickle.dump(archive_info, pklh)
            exit(1)
        os.chdir(dir_run)
    else:
        logging.info('Exited without archiving')
        exit(0)
        
@click.command()
@main_opt_get('d')
@main_opt_get('k')
@main_opt_get('K')
@main_opt_get('t')
@main_opt_get('u')
@main_opt_get('f')
def main(delete, keep, keep_config, keep_tarball, uncompressed, files):
    archive(delete, keep, keep_config, keep_tarball, compression=(not uncompressed), files=files)

@click.command()
@main_opt_get('d')
@main_opt_get('k')
@main_opt_get('K')
@main_opt_get('t')
@main_opt_get('u')
@main_opt_get('f')
def safe(delete, keep, keep_config, keep_tarball, uncompressed, files):
    archive(delete, keep, keep_config, keep_tarball, True, compression=(not uncompressed), files=files)

if __name__ == '__main__':
    main()
