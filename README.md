# Archiving script for the Goate Lab

This script does the following:

*  Checks if the directory is ready to archive:
   *  Makes sure you are using screen or tmux so that the archiving doesn't fail if you get disconnected
   *  Checks if there are uncommitted or local changes to any git repositories
*  Compresses all of the files in the current directory into a tarball.
   *  Uses a parallelized cluster job for speed
   *  Uses the user's scratch directory to conserve space
   *  Checks that the archiving was successful and reruns if necessary
*  *Optionally* deletes either all files or all large files after creating the tarball to save space on the server.
*  Archives the tarball to tape using the dsmc command
*  Deletes the tarball (optional)
*  Saves the list of archived files along with which ones were deleted, sizes, and permissions to a TSV file in the directory and to a lab database
   *  Also records who did the archiving, the date and time of archiving and the current path to a log file and the database
*  Saves an extensive log file of the archiving operations
   
## Installing the script:

Install the script and dependancies using the following command:

```bash
mamba install -c bfh -c conda-forge archive_goate
```

## Running the script:

First, make sure you are using `screen` or `tmux` so that the archiving doesn't fail if you get disconnected.

***NOTE:*** `dsmc` only works on the login nodes, so make sure you are not in an `lsf` job when using the script. That means no `ijob`, no `regjob`, no `bsub`, etc.

1.  To start a new `screen` session, run `screen` from the command line
2.  To detach from your session, press `Ctrl+a` then press `d`
3.  To resume your session if you are disconnected or detach, run `screen -r` from the command line

You can type `archive_goate` in the directory you want to archive. If you definitely don't want to delete files and you want a safe script, run `archive_goate_safe`.

The scripts now take arguments for what you want to do during the archiving process:

```
Usage: archive_goate [OPTIONS]

Options:
  -d, --delete                    Delete files after confirming
  -k, --keep [small|none|ask|default]
                                  Keep small files, none, or ask (default is
                                  ask when deleting)
  -K, --keep-config FILENAME      Use YAML file to determine which files to
                                  keep
  -t, --keep-tarball [yes|no|ask]
                                  Keep tarball after archiving
  -u, --uncompressed              Do not compress tarball with bzip2
  -f, --files TEXT                Files to archive
  --help                          Show this message and exit.
```

The script will prompt you to determine how to proceed:

1.  The script will normally scan your current directory for files and archive everything. If you have large files in the directory and want to archive specific files, you can use `-f` or `--files`, one time for each file.
1.  If you are not using screen or tmux, it will ask you if you want to quit and use either one.
1.  If any unsaved changes exist in git tracked directories and have specified `-d` or `--delete` on the command line, it will show the status and ask if you want to quit to commit/push changes
1.  It will show the archiving info and ask if you want to proceed with the process.
1.  The script will generally compress files when making the tarball. If the files you are archiving are incompressible (already compressed), you can select `-u` or `--uncompressed` when calling the script.
1.  After generating the tarball, it will confirm if you want to delete files if you have specified `-d` or `--delete` on the command line.
    *  If you say yes, it will ask if you want to keep small files unless you have specified using `-k`/`--keep` or `-K`/`--keep-config`.
1.  Before archiving, it will ask if you want to keep the tarball unless specified using `-t` or `--keep-tarball`. This is not recommended unless deleting files because it can as much as double space usage.

If there are errors during compression or archiving, the script will detect them and ask if you want to proceed, try again or cancel.

### Small files

You can choose to keep small files if you are deleting files while you archive. By default, the following files are kept:

*  Logs ('.log', '.err', '.stderr', '.out', '.stdout') under 200 MiB
*  Config files ('.yaml', '.yml', '.conf', '.json', '.ini') under 50 MiB
*  Scripts ('.sh', '.r', '.rmd', '.smk', '.py') under 5 MiB
*  Text and output files ('.txt', '.pdf', '.html', '.htm', '.org', '.md') under 100 MiB
*  Images ('.tiff', '.png', '.jpg', '.svg') under 5 MiB
*  Readme files ('.readme') under 10 MiB

You can override those settings by using a yaml file specified like `--keep-config archive.yaml` in the directory with the following format:

```yaml
types:
  Type1:
    - '.ext1'
    - '.ext2'
  Type2:
    - '.ext3'
    - '.ext4'
    - '.ext5'
sizes:
  Type1: 100
  Type2: 10
```

Sizes are in MiB.

### Specified files

Sometimes, you will only want to archive specific files. You can select them with `-f` or `--files`, one time for each file. So to archive `file1` and `file2`, do the following:

```bash
archive_goate -f file1 -f file2
```

This is particularly useful in batch mode.

### Batch mode

When using `archive_goate_safe`, batch mode is available. This will allow you to use the script with no prompts, such as in a loop.

***MAKE SURE YOU ARE USING TMUX OR SCREEN WHEN USING THIS MODE BECAUSE THE SCRIPT WILL NOT CHECK***

An example of how to use the script for two already-compressed files follows:

```bash
archive_goate_safe --batch --uncompressed --keep-tarball no -f file1 -f file2
```

It becomes much more useful to use this in a loop in the shell of your choice. The following code archives sets of cram and crai files:

```bash
# Only run the following line the first time to initialize the log with completed files
echo start_of_file > archiving.done.log

# This is the main loop
while read f; do # lines in ls get assigned to $f
  success=0 # set success to good exit code
  if grep -qv $f archiving.done.log; then # check if file already archived
    archive_goate_safe -f $f -f $f.crai -t no -u --batch # main command
    success=$? # store exit code for script
  fi
  if [ $success -ne 0 ]; then
    break # exit loop if unsuccessful
  else
    echo $f >> archiving.done.log # add to completion log if successful
  fi
done < <(ls -1 *.cram) # search for cram files
```

Once you have finished running this sample script, you can use another loop or the `xargs` command to delete all files in `archiving.done.log`.

## Importing the functions in Python or snakemake

To use the archiving script in snakemake, the `archive_snakemake.py` script is provided in this repo. You must use snakemake with the local profile or as a localrule (not on the lsf cluster), and `archive_goate` must be in your environment. An archiving rule looks like this:

```
rule archive_files:
    input:
        archive = [list of files to be archived],
        prereqs = [list of files to wait for]
    params:
        directory = [directory to use as root for the archiving if
            that differs from the common prefix of the input files],
        delete = [if True, the input files will be deleted after archiving],
        keep = [if 'small', small files will be kept after archiving; optional],
        compress = [if False, the tarball will not be compressed; optional],
    localrule: true
    output: [output file to be touched after archiving]
    script: 'archive_snakemake.py'
```

To directly use the script in python, just import it:

```python
import os
import archive_goate

DIRECTORY = "/path/to/files"
FILES = ["file1", "file2"]

DELETE = False
KEEP = "none"
SAFE = True
COMPRESS = True

directory_norm = os.path.normpath(DIRECTORY)

archive_goate.archive(DELETE, KEEP, safe=SAFE, batch=True,
                      files=FILES, directory=directory_norm,
                      compression=COMPRESS)
```
