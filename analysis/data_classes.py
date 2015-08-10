from __future__ import division
__author__ = 'chris'


import tables as tb
import odor_analysis
from sniff import add_sniff_events, SniffExistExemption


class Recording(tb.File):
    """
    Container for single recording worth of data. Inherits tables file object, and provides container for methods that
    can be used to work with data within these files.
    """


    # ----  DECLARE CLASS FUNCTIONS. -----
    add_sniff_events = add_sniff_events

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
    plot_odor_rasters_sniff = odor_analysis.plot_odor_rasters_sniff

    def __init__(self, filename, *args, **kwargs):

        super(Recording, self).__init__(filename, mode='r+', *args, **kwargs)  # open file in append mode.
        assert hasattr(self.root, 'clusters')
        assert hasattr(self.root, 'events')
        assert hasattr(self.root, 'streams')
        try:
            add_sniff_events(self, overwrite=False)
        except SniffExistExemption:  # if sniff events already exists, don't worry about it!
            pass

        return

    def get_stream(self, stream, start, end, timebase='event'):
        """
        This returns a stream in the event timebase (ie. for plotting).

        :param stream: a the stream that want to work with. This can be a streamname or a stream node.
        :param start: the start sample that you want returned. This is in the event timebase.
        :param end: the end sample that you want returned. This is in the event timebase

        :returns: stream[start:end]
        """

        if isinstance(stream, str):
            stream = self.get_node('/streams/', stream)

        ev_timebase = self.root.events._v_attrs['sample_rate_Hz']
        st_timebase = stream._v_attrs['sample_rate_Hz']
        conversion_factor = ev_timebase / st_timebase

        start_streambase = start / conversion_factor
        end_streambase = end / conversion_factor

        return stream[start_streambase:end_streambase]

    def sample_to_ms(self, samples):
        """
        Convert sample time value(s) to millisecond timebase.

        :param samples: numeric (array or scalar)
        """

        sample_rate = self.root.events._v_attrs['sample_rate_Hz']
        ms = samples/sample_rate * 1000
        return ms


class Recordings(list):
    """
    Container for multiple recordings. Functions here override the original functions to work on the whole list of
    recordings with a single function call.
    """
    pass
