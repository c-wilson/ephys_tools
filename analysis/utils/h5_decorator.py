__author__ = 'chris'
import tables as tb


class h5decorator(object):
    """
    This little guy checks whether an h5 parameter is a string. If it is, it runs the function by opening the h5.
    """

    def __init__(self, f):
        self.f = f
        return

    def __call__(self, h5, *args, **kwargs):
        if isinstance(h5, str):
            with tb.open_file(h5) as h5:
                return self.f(h5, *args, **kwargs)
        else:
            return self.f(*args, **kwargs)