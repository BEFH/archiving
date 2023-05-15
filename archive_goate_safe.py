#!/usr/bin/env python3

__version__ = '2.3'

import logging
import pickle
import os

try:
    import archive_goate as arch
except:
    import archive_minerva as arch

def main():
    logging.basicConfig(
        filename='archive_{}.log'.format(arch.logtime),
        level=logging.INFO,
        format='%(asctime)s %(levelname)s: %(message)s',
        datefmt='%m/%d/%Y %H:%M:%S')

    arch.check_tmux()

    settings = arch.load_config('archive.yaml')

    print('Scanning for files')
    assert(__version__ == arch.__version__)
    logging.info(f'Archiving script v{__version__}')
    logging.info('Scanning directory for files')
    files = arch.list_files(settings)
    logging.info('Done scanning directory for files')

    if any(files['kind'] == 'gitdir'):
        for gitdir in files[files['kind'] == 'gitdir']['path']:
            arch.git_test_and_prompt(gitdir)

    log_sz, sizes = arch.get_sizes(files)

    if arch.question('Archive this folder:\n{}?'.format(os.getcwd())):
        archive_info, temp_tarball = arch.make_tarball(files, sizes['total'])
        files, archive_info = arch.delete_no_files(files, archive_info, sizes)
        arch.tsm_archive(temp_tarball, archive_info, files)
        database = '/sc/arion/projects/LOAD/archive/archive.sqlite'
        arch.write_database(database, files, archive_info)
        arch.write_tables(files, archive_info)
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
