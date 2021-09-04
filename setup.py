from setuptools import setup, find_packages

setup(
    name='archive_goate',
    version='1.0b4',
    description='Automatic backup script for goate lab.',
    long_description=('Using DSMC to back up files is somewhat compled and\n'
                      'it\'s difficult to track what files and directories\n'
                      'have been archived in the lab. This script automates\n'
                      'the process of archiving, stores information about the\n'
                      'archive, and stores info across users for the entire\n'
                      'lab.'),
    url='https://github.com/BEFH/archiving',
    author='Brian Fulton-Howard',
    author_email='brian.fulton-howard@mssm.edu',
    license='MIT',
    keywords='archive Linux TSM dsmc tar',
    py_modules=['archive_goate'],
    classifiers=[
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
        'Operating System :: POSIX :: Linux',
    ],
    install_requires=['pytz','pandas'],
    entry_points={
        'console_scripts': ['archive_goate=archive_goate:main'],
    },
)
