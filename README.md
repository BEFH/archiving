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
*  Deletes the tarball
*  Saves the list of archived files along with which ones were deleted, sizes, and permissions to a TSV file in the directory and to a lab database
   *  Also records who did the archiving, the date and time of archiving and the current path to a log file and the database
*  Saves an extensive log file of the archiving operations
   
## Installing the script:

Install the script and dependancies using the following command:

```bash
mamba install -c befh -c conda-forge archive_goate
```

## Running the script:

First, make sure you are using `screen` or `tmux` so that the archiving doesn't fail if you get disconnected.

1.  To start a new `screen` session, run `screen` from the command line
2.  To detach from your session, press `Ctrl+a` then press `d`
3.  To resume your session if you are disconnected or detach, run `screen -r` from the command line

You can type `archive_goate` in the directory you want to archive. If you definitely don't want to delete files and you want a safe script, run `archive_goate_safe`.

The script will prompt you to determine how to proceed:

1.  If you are not using screen or tmux, it will ask you if you want to quit and use either one.
1.  If any unsaved changes exist in git tracked directories, it will show the status and ask if you want to quit to commit/push changes
1.  It will show the archiving info and ask if you want to proceed with the process.
1.  After generating the tarball, it will ask if you want to delete files.
    *  If you say yes, it will ask if you want to keep small files

If there are errors during compression or archiving, the script will detect them and ask if you want to proceed, try again or cancel.

### Small files

You can choose to keep small files if you are deleting files while you archive. By default, the following files are kept:

*  Logs ('.log', '.err', '.stderr', '.out', '.stdout') under 200 MiB
*  Config files ('.yaml', '.yml', '.conf', '.json', '.ini') under 50 MiB
*  Scripts ('.sh', '.r', '.rmd', '.smk', '.py') under 5 MiB
*  Text and output files ('.txt', '.pdf', '.html', '.htm', '.org', '.md') under 100 MiB
*  Images ('.tiff', '.png', '.jpg', '.svg') under 5 MiB
*  Readme files ('.readme') under 10 MiB

You can override those settings by putting a file called `archive.yaml` in the directory with the following format:

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

