#!/usr/bin/env python3
'''
Snakemake script to archive files using the archive_goate module.

Usage:
    input:
        archive = [list of files to be archived],
        prereqs = [list of files to wait for]
    params:
        directory = [directory to use as root for the archiving if
            that differs from the common prefix of the input files],
        delete = [if True, the input files will be deleted after archiving],
        keep = [if 'small', small files will be kept after archiving; optional],
    output: [output file to be touched after archiving]
'''

import os
import re
import archive_goate


if 'snakemake' not in locals():
    class snakemake:  # pylint: disable=C0103, R0903
        '''
        This is a mock class to be used in the snakemake script.
        '''
        input = {'archive': ['/sc/arion/projects/LOAD/archive/admin/sandbox/' + x
                             for x in ['archive_10-Nov-2023_13.29.log',
                                       'archive_test_1/archive_dump.p']]}
        output = ['archive_done.txt']
        params = {'delete': True,
                  'directory': '/sc/arion/projects/LOAD/archive/'}


# Get the input files (unique inputs)
files = list(set(snakemake.input['archive']))

if "directory" in snakemake.params:
    directory = snakemake.params["directory"]
    directory_norm = os.path.normpath(directory)
else:
    directory_norm = "."  # pylint: disable=C0103

if "delete" in snakemake.params:
    DELETE = snakemake.params["delete"]
    SAFE = False
else:
    DELETE = False
    SAFE = True

if "keep" in snakemake.params:
    KEEP = snakemake.params["keep"]
    if DELETE and KEEP not in ("small", "none"):
        raise ValueError("params['keep'] must be 'small' or 'none' if "
                         "params['delete'] is True")
    if not DELETE and KEEP != "default":
        raise ValueError("params['keep'] must be 'default' if "
                         "params['delete'] is False")
elif DELETE:
    KEEP = "none"
else:
    KEEP = "default"

# Get common prefix
prefix = os.path.commonprefix([os.path.dirname(x) for x in files])
prefix_norm = os.path.normpath(prefix)

for f in files:
    if not os.path.exists(f):
        raise FileNotFoundError(f)

if all(f.startswith("/") for f in files):
    # Remove prefix from files
    files = [re.sub("^" + prefix_norm + "/", "", x) for x in files]
    if "directory" in locals() and directory_norm != prefix_norm:
        if not re.search("^" + directory_norm, prefix_norm):
            raise ValueError("params['directory'] not in input paths")
        extra_path = re.sub("^" + directory_norm + "/", "", prefix_norm)
        files = [os.path.join(extra_path, x) for x in files]
    else:
        directory_norm = prefix_norm
elif directory in locals() and directory_norm != '.':
    if not re.search("^" + directory_norm, prefix_norm):
        raise ValueError("params['directory'] not in input paths")
    files = [re.sub("^" + directory_norm + "/", "", x) for x in files]

archive_goate.archive(DELETE, KEEP, safe=SAFE, batch=True,
                      files=files, directory=directory_norm)

# Touch output file
with open(snakemake.output[0], 'w', encoding='utf-8') as f:
    print("Archiving done.", file=f)
