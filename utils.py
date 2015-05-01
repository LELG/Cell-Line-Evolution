"""
Simple utility functions which are used
throughout the simulation.

Authors
-------
Yoshua Wakeham : yoshua.wakeham@petermac.org


Date Created
------------
4 Feb 2015


Notes
-----
The module is designed to be used with Python 2.7.x.


Change Log
----------
4 Feb 2015: created secs_to_hms() utility
            moved directory creation utility here
1 May 2015: created delete_local_file() utility
"""
import os
import errno


def secs_to_hms(time_delta):
    """Convert a time delta in seconds to HH:MM:SS.S representation."""
    mins, secs = divmod(time_delta, 60)
    hours, mins = divmod(mins, 60)
    return "{:02d}:{:02d}:{:04.1f}".format(int(hours), int(mins), secs)


def make_path_unless_exists(path):
    """Make a directory, unless it already exists"""
    try:
        os.makedirs(path)
    except OSError as exception:
        if exception.errno != errno.EEXIST:
            raise


def delete_local_file(fname):
    """Delete a file from the local directory (if it exists).

    This try-except is taken from the StackOverflow post:
        'Most pythonic way to delete a file which may not exist'
    """
    try:
        os.remove(fname)
    except OSError as err:
        if err.errno != errno.ENOENT:
            raise
