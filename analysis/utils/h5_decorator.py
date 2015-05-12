__author__ = 'chris'
import tables as tb


def h5decorator(f):
    """
    This little guy checks whether an h5 parameter is a string. If it is, it runs the function by opening the h5.

    The big benefit here is that it allows general functions that can accept both open files and filenames. Since tables
    objects don't pickle, passing filenames allows for files to be opened and processed when using parallel/cluster
    processing.
    """

    def deco(h5, *args, **kwargs):
        if isinstance(h5, str):
            with tb.open_file(h5) as h5:
                return f(h5, *args, **kwargs)
        else:
            return f(h5, *args, **kwargs)
    return deco