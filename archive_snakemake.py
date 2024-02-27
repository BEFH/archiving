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
        compress = [if False, the tarball will not be compressed; optional],
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
        input = {'archive': ['/sc/arion/projects/load/data-ext/ADSP/ADSP_Umbrella_ng00067.v10/CRAMs/snowball/fromsnowball/distribution/adsp/cram/snd10001/' + x
                             for x in ['DS-AGEADLT-IRB-PUB/A-ACT-AC003412-BR-NCR-16AD84907_vcpa1.0_DS-AGEADLT-IRB-PUB.cram',
                                       'DS-AGEADLT-IRB-PUB/A-ACT-AC003412-BR-NCR-16AD84907_vcpa1.0_DS-AGEADLT-IRB-PUB.cram.crai']]}
        output = ['archive_done.txt']
        params = {'delete': False,
                  'directory': '/sc/arion/projects/load/data-ext/ADSP/ADSP_Umbrella_ng00067.v10/CRAMs/snowball/fromsnowball/distribution/',
                  'compress': True}


# Get the input files (unique inputs)
files = list(set(snakemake.input['archive']))

params = dict(snakemake.params)

if "directory" in params:
    directory = params["directory"]
    directory_norm = os.path.normpath(directory)
else:
    directory_norm = "."  # pylint: disable=C0103

if "delete" in params:
    DELETE = params["delete"]
    SAFE = False
else:
    DELETE = False
    SAFE = True

if "keep" in params:
    KEEP = params["keep"]
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

if "compress" in params:
    COMPRESS = params["compress"]
else:
    COMPRESS = True

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

#import ipdb; ipdb.set_trace()

archive_goate.archive(DELETE, KEEP, safe=SAFE, batch=True,
                      files=files, directory=directory_norm,
                      compression=COMPRESS)

# Touch output file
with open(snakemake.output[0], 'w', encoding='utf-8') as f:
    print("Archiving done.", file=f)
