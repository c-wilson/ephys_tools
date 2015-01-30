from __future__ import division
__author__ = 'chris'


import tables as tb
import numpy as np
import odor_analysis
from sniff import make_sniff_events, SniffExistExemption


class Recording(tb.File):
    """
    Container for single recording worth of data. Inherits tables file object, and provides container for methods that
    can be used to work with data within these files.
    """


    # ----  DECLARE CLASS FUNCTIONS. -----
    make_sniff_events = make_sniff_events

    # ---- event finders
    find_odor_events = odor_analysis.find_odor_events
    find_sniff_events = odor_analysis.find_sniff_events

    # ---- Odor functions
    get_unique_odors = odor_analysis.get_unique_odors
    get_unique_concentrations = odor_analysis.get_unique_concentrations
    get_odor_name = odor_analysis.get_odor_name

    # ---- PSTH/Raster functions
    get_odor_psth = odor_analysis.get_odor_psth
    get_rasters = odor_analysis.get_rasters

    # ---- PLOTTING Functions
    plot_odor_psth_no_baseline = odor_analysis.plot_odor_psth_no_baseline
    plot_odor_psth_w_baseline = odor_analysis.plot_odor_psth_w_baseline
    plot_odor_rasters = odor_analysis.plot_odor_rasters

    def __init__(self, filename, *args, **kwargs):

        super(Recording, self).__init__(filename, mode='a', *args, **kwargs)  # open file in append mode.
        assert hasattr(self.root, 'clusters')
        assert hasattr(self.root, 'events')
        assert hasattr(self.root, 'streams')
        try:
            make_sniff_events(self, overwrite=False)
        except SniffExistExemption:  # if sniff events already exists, don't worry about it!
            pass

        return

# TODO: make recordings class.
class Recordings(list):
    """
    Container for multiple recordings. Functions here override the original functions to work on the whole list of
    recordings with a single function call.
    """
    pass