#!/usr/bin/env python
# Archiving script Version 1.1.2

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
from itertools import compress

# packages
import pytz
import pandas as pd


logtime = datetime.datetime.now().strftime('%d-%b-%Y_%H.%M')

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
    path = entry.path
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
    path = entry.path
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


def git_test_and_prompt(gitdir):
    git_status = git_test(gitdir)
    if git_status == 1:
        print("We recommend you cancel the script and commit changes.")
        if question("Do you want to cancel?", True):
            exit(1)
        else:
            logging.warning(
                '{} is a git project and is ahead of remote'.format(gitdir))
    elif git_status == 2:
        if question("Do you want to cancel to commit and push changes?", True):
            exit(1)
        else:
            logging.warning(
                '{} is a git project and has untracked changes'.format(gitdir))
    else:
        logging.info('{} is a git project and is up to date'.format(gitdir))


def question(text, default=None, prepend=''):
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
        question(text, default, 'Please enter a response:\n')
    elif answer.lower()[0] in ["y", "yes"]:
        return True
    elif answer.lower()[0] in ["n", "no"]:
        return False
    else:
        question(text, default, 'Invalid answer; please try again:\n')


def scantree(path, settings):
    '''Recursively yield DirEntry objects for given directory.'''
    for entry in os.scandir(path):
        if entry.is_dir(follow_symlinks=False):
            if entry.name == '.git':
                yield getinfo(entry, settings)
            yield from scantree(entry.path, settings)
        else:
            yield getinfo(entry, settings)


def list_files(settings):
    files = [x for x in scantree('.', settings) if x is not None]
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
    return(files)


def check_tmux():
    '''Check that screen or tmux are in use:'''
    if os.getenv('TERM') != 'screen':
        print('We recommend running in screen or tmux.')
        if question('Quit to run screen or tmux?', True):
            exit(2)


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


def get_sizes(files):
    files_nobroken = files[files['mode'].notnull()]
    files_nobroken = files_nobroken.astype({'mode': 'int32'})

    total = files.loc[:, 'size_mib'].sum().round(1)
    freed = files[files['keep'] == 'no'].loc[:, 'size_mib'].sum().round(1)
    kept = (total - freed).round(1)

    kinds = files_nobroken.groupby(['kind']).sum().round(3)[['size_mib']]

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


def make_tarball(files, total):
    fullpath = os.path.normpath(os.getcwd())
    parent = os.path.dirname(fullpath)
    cwd = os.path.basename(fullpath)
    tznm = '/'.join(os.readlink('/etc/localtime').split('/')[-2:])
    dt_current = datetime.datetime.now(tz=pytz.timezone(tznm))
    date = dt_current.strftime('%G-%b-%d')
    dt_iso = dt_current.isoformat()
    hash_dtd = filehash(cwd, dt_iso)
    files['file_id'] = [filehash(x, dt_iso) for x in files['path']]
    files['archive_id'] = hash_dtd
    tarball = '{dirname}_{date}.tar.bz2'.format(dirname=cwd, date=date)
    username = getpass.getuser()

    archive_info = {
        'archive_id': hash_dtd,
        'user_name': username,
        'time': dt_current,
        'archive_name': tarball,
        'archive_directory': fullpath,
        'total_mib': total
        }

    oo = '{dirname}_{date}.tarball.stdout'.format(dirname=cwd, date=date)
    eo = '{dirname}_{date}.tarball.stderr'.format(dirname=cwd, date=date)

    tempdir = '/sc/arion/scratch/{}/archiving'.format(username)
    temp_tarball = '{}/{}'.format(tempdir, tarball)
    tar_cmd = 'tar -cvpjf {} -C {} {}'.format(temp_tarball, parent, cwd)
    tar_cmd_fast = 'tar c -I"pbzip2 -p24" -f {} -C {} --xattrs --acls {}'
    tar_cmd_fast = tar_cmd_fast.format(temp_tarball, parent, cwd)

    try:
        not_used = subprocess.run(["pbzip2", "-h"], capture_output=True)
    except FileNotFoundError:
        pbgzip = False
    else:
        pbzip = True
        
    newtar = test_tar()

    ncore = 1
    if newtar and pbzip:
        logging.info("Using pbzip2 to compress in parallel. ACLs and XATTRs recorded.")
        tar_cmd = tar_cmd_fast
        ncore = 24
    elif pbzip:
        print("Tar is older than v1.28, so multithreading and ACLs are unsupported.")
        print("Use 'mamba install -c conda-forge tar' to enable features.")
        logging.warning("Tar is older than 1.28. Using serial compression and no ACLs")
    elif newtar:
        print("Parallel bzip2 missing. Archiving will be slow unless you install")
        print("pbzip2. Use 'mamba install -c conda-forge pbzip2' to enable.")
        logging.warning("pbzip2 is missing. Using slow serial compression.")
    else:
        print("Parallel bzip2 missing and tar is older than v1.28 so multithreading")
        print("and ACLs are unsupported.")
        print("Use 'mamba install -c conda-forge pbzip2 tar' to enable features.")
        logging.warning("pbzip2 is missing and tar is older than 1.28.")
        logging.warning("Using slow serial compression and no ACLs.")

    job_cmd = ('bsub -P acc_LOAD -q premium -n {} -R rusage[mem=1000] '
               + '''-R span[hosts=1] -W 140:00 -oo {} -eo {} '{}' ''')
    job_cmd = job_cmd.format(ncore, oo, eo, tar_cmd)
    # -c create
    # -v verbose
    # -p preserve permissions
    # -j bzip2
    # -f tarball file name
    tempdir_path = pathlib.Path(tempdir)
    pathlib.Path.mkdir(tempdir_path, parents=True, exist_ok=True)
    archive_job = subprocess.run(job_cmd, capture_output=True, shell=True)
    print(archive_job.stdout.decode())
    print(archive_job.stderr.decode())
    archive_job.check_returncode()
    job_expression = r'(?<=^Job <)\d+(?=> is submitted to queue)'
    job = re.search(job_expression, archive_job.stdout.decode())[0]
    completed = False
    while not completed:
        time.sleep(10)
        jobstat = subprocess.run('bjobs -do "stat" -noheader {}'.format(job),
                                 capture_output=True, shell=True)
        completed = jobstat.stdout.decode().strip() not in ['RUN', 'PEND']
    assert jobstat != 'EXIT', "ERROR: Archiving failed!"
    return archive_info, temp_tarball


def filehash(file, date):
    '''Generate hash from file/dir name and date'''
    hash = hashlib.shake_256()
    hash.update((file + date).encode())
    return hash.hexdigest(8)


def file_rm(files, archive_info, sizes, log_sz):
    logging.info('Removing files')
    files['removal'] = 'not removed'
    if question('Do you want to keep {} MiB of small files?'.format(
                sizes['kept']), True):
        logging.info('Keeping {} MiB of small files'.format(sizes['kept']))
        logprint(log_sz)
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


def delete_files(files, archive_info, sizes, log_sz):
    idx = files[files['filename'] == 'archive_{}.log'.format(logtime)].index
    files.drop(idx, inplace=True)
    if question('Do you want to remove files after generating archive', False):
        try:
            files, archive_info = file_rm(files, archive_info, sizes, log_sz)
        except:
            logging.exception('Uncaught exception deleting files')
            print("Unexpected error deleting files:", sys.exc_info()[0])
            with open('archive_dump.p', 'wb') as pklh:
                pickle.dump(files, pklh)
            if question('Do you want to continue anyway?', False):
                files['removal'] = 'fail: other'
                archive_info['freed_mib'] = 0
                archive_info['kept_mib'] = sizes['total']
                archive_info['removal_failure'] = 'yes'
            else:
                archive_info['removal_failure'] = 'yes'
                with open('info_dump.p', 'wb') as pklh:
                    pickle.dump(archive_info, pklh)
                raise
    else:
        logging.info('Not removing files')
        archive_info['freed_mib'] = 0
        archive_info['kept_mib'] = sizes['total']
        files['removal'] = 'kept'
    removed = files.loc[files['removal'] == 'success', 'size_mib'].sum()
    removed = removed.round()
    discrep = archive_info['freed_mib'].round() - removed
    removed = 'Removed {} MiB of files.'.format(int(removed))
    print(removed)
    logging.info(removed)
    if discrep != 0:
        discrep = 'Could not remove {} MiB of files.'.format(int(discrep))
        print(discrep)
        logging.warning(discrep)
    return files, archive_info


def tsm_archive(temp_tarball, archive_info, files=None, attempt=1):
    archive_name = archive_info['archive_name']
    shutil.copy(temp_tarball, archive_name)
    os.remove(temp_tarball)
    logging.info('DSMC Job starting')
    print('DSMC Job starting')
    dsmc_cmd = 'dsmc archive -se={} "{}"'.format(archive_info['user_name'], archive_name)
    dsmc_job = subprocess.run(dsmc_cmd, capture_output=True, shell=True)
    print(dsmc_job.stdout.decode())
    print(dsmc_job.stderr.decode())
    logging.info('DSMC Job completed')
    logprint(dsmc_job.stdout.decode())
    logprint(dsmc_job.stderr.decode())
    try:
        dsmc_job.check_returncode()
    except:
        if attempt < 4 and question('Archiving failed with {}. Try again?'.format, False):
            tsm_archive(temp_tarball, archive_info, files, attempt + 1)
        else:
            if files is not None:
                with open('archive_dump.p', 'wb') as pklh:
                    pickle.dump(files, pklh)
            archive_info['exception'] = 'archiving'
            with open('info_dump.p', 'wb') as pklh:
                pickle.dump(archive_info, pklh)
            logging.exception('Uncaught exception archiving files')
            logging.exception('Please contact Brian to continue the archiving.')
            print("Unexpected error archiving files:", sys.exc_info()[0])
            print("\033[1m\033[91mContact Brian to finish archiving!\033[0m")
            raise
    else:
        logging.info('DSMC Job successful')
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


def write_tables(files, archive_info):
    try:
        fname = 'archived_file-list_{}.tsv.gz'.format(logtime)
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
        logprint(archive_info)
        with open('archive_info_{}.log'.format(logtime), 'w') as ai:
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


def logprint(message):
    with open('archive_{}.log'.format(logtime), 'a') as logfile:
        print('\n{}\n'.format(message), file=logfile)


def main():
    logging.basicConfig(
        filename='archive_{}.log'.format(logtime), level=logging.INFO,
        format='%(asctime)s %(levelname)s: %(message)s',
        datefmt='%m/%d/%Y %H:%M:%S')

    check_tmux()

    settings = load_config('archive.yaml')

    print('Scanning for files')
    logging.info('Archiving script v1.1.2')
    logging.info('Scanning directory for files')
    files = list_files(settings)
    logging.info('Done scanning directory for files')

    if any(files['kind'] == 'gitdir'):
        for gitdir in files[files['kind'] == 'gitdir']['path']:
            git_test_and_prompt(gitdir)

    log_sz, sizes = get_sizes(files)

    if question('Archive this folder:\n{}?'.format(os.getcwd())):
        archive_info, temp_tarball = make_tarball(files, sizes['total'])
        files, archive_info = delete_files(files, archive_info, sizes, log_sz)
        tsm_archive(temp_tarball, archive_info, files)
        database = '/sc/arion/projects/LOAD/archive/archive.sqlite'
        write_database(database, files, archive_info)
        write_tables(files, archive_info)
        logging.info('Done')
        if 'removal_failure' in archive_info:
            with open('info_dump.p', 'wb') as pklh:
                pickle.dump(archive_info, pklh)
            exit(1)
        exit(0)
    else:
        logging.info('Exited without archiving')
        exit(0)


if __name__ == '__main__':
    main()
