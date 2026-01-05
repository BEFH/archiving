#!/usr/bin/env python3

"""
This script is an admin tool to prune backups of the SQL database.
"""

import os
import time
from datetime import datetime, timedelta
import click

ARCHIVE_DIR = '/sc/arion/projects/LOAD/archive/admin/bk'
NOW = time.time()

def mtime(path):
    path = os.path.join(ARCHIVE_DIR, path)
    return os.stat(path).st_mtime

def is_friday(ts):
    return datetime.fromtimestamp(ts).weekday() == 4

def is_last_friday(ts):
    d = datetime.fromtimestamp(ts)
    return d.weekday() == 4 and (d + timedelta(days=7)).month != d.month

@click.command()
def main():
    files = [f for f in os.listdir(ARCHIVE_DIR) if f.endswith('.sql.gz')]
    files.sort(key=mtime, reverse=True)

    keep = {}

    for f in files[-2:]:
        d = datetime.fromtimestamp(mtime(f))
        keep[('oldest', d.strftime('%Y-%m-%d'))] = f

    for f in files[:-2]:
        ts = mtime(f)
        age_days = int((NOW - ts) / 86400)
        d = datetime.fromtimestamp(ts)

        key = None
        if age_days <= 3:
            keep[('all', d.strftime('%Y-%m-%d_%H:%M'))] = f
            continue
        elif age_days <= 14:
            key = ('daily', d.strftime('%Y-%m-%d'))
        elif age_days <= 48 and is_friday(ts):
            key = ('weekly', d.strftime('%G-%V'))
        elif age_days > 48 and is_last_friday(ts):
            key = ('monthly', d.strftime('%Y-%m'))

        if key and key not in keep:
            keep[key] = f

    to_delete_file = [f for f in files if f not in keep.values()]

    print("=== WILL KEEP ===")
    for k, f in keep.items():
        reason, date = k
        print(f'keeping {reason} from {date}: {f}')

    if not to_delete_file:
        print("\nNo files to delete.")
        return (keep, [])
    else:
        print("\n=== WILL DELETE ===")
        for f in to_delete_file:
            print(f)

    if click.confirm("\nDo you want to delete these files?", default=False):
        to_delete_path = [os.path.join(ARCHIVE_DIR, f) for f in to_delete_file]
        for f in to_delete_path:
            os.remove(f)
        print(f"Deleted {len(to_delete_path)} files.")
    else:
        print("Aborted.")

    return (keep, to_delete_file)

if __name__ == '__main__':
    keep, delete = main()