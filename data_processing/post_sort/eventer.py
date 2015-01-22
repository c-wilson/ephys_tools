from __future__ import division
__author__ = 'chris'

"""
This is a module for making event structures consistent between Voyeur and output from the electophysiology system
"""

import numpy as np
import tables as tb
from serial_parser import parse_serial_stream
import logging
try:
    import cPickle as pickle  # faster, better.
except ImportError:
    import pickle

def eventer(raw_kwd_fn, destination):
    """

    :param raw_kwd_fn:
    :param destination:
    :return:
    """

    assert isinstance(destination, tb.File)
    with tb.open_file(raw_kwd_fn, 'r') as raw_kwd:
        logging.info('Processing trial events.')
        make_trial_upload_events(raw_kwd, destination)
        logging.info('Processing run events.')
        make_run_events(raw_kwd, destination)
        for handler in handlers:
            logging.info(u'Making events for {0:s}'.format(handler.__name__))
            handler(raw_kwd, destination)
    logging.info('Event building complete.')


def make_trial_upload_events(raw_kwd, dest_file):
    """

    :param raw_kwd: hdf5 object with raw data from process.py
    :param dest_file: hdf5 object to contain processed data.
    :return:
    """

    assert isinstance(raw_kwd, tb.File)
    assert isinstance(dest_file, tb.File)
    grp = dest_file.create_group(u'/events',
                                 u'trials',
                                 title=u'Trial parameters uploaded to arduino from Voyeur.',
                                 createparents=True)
    n_events = 0
    trial_sets = []
    serial_dicts = []
    fieldsets = []
    sample_offsets = []
    offset = 0

    for r in xrange(raw_kwd.root.recordings._v_nchildren):
        rec = raw_kwd.get_node(u'/recordings/{0:d}'.format(r))
        serial_st = rec.serial_trial
        try:
            fs = serial_st.get_attr('sample_rate_Hz')
        except AttributeError:
            fs = 20833.
            logging.warning(u'Warning, no sampling rate specified for serial trial number. Default is set to 20833 Hz.')
        serial_dict = parse_serial_stream(serial_st, fs=fs, word_len=2, baudrate=300)

        n_events += len(serial_dict)
        serial_dicts.append(serial_dict)

        # do some parity check to make sure that there are no repeated trial numbers and that the trial numbers found in
        # the Voyeur file are also found in the serial stream.
        v_trials = rec.Voyeur_data.Trials
        v_trial_nums = v_trials[:]['trialNumber']
        s_trial_nums = serial_dict.keys()
        if len(v_trial_nums) != len(np.unique(v_trial_nums)):
            # trial number are duplicated within the voyeur table, find which trials are duplicates and raise exemption.
            t, n = np.unique(v_trial_nums, return_counts=True)
            repeats = t[n > 1]
            raise ValueError('Voyeur table has repeated trial numbers: {0}.'.format(repeats))
        if len(s_trial_nums) != len(np.unique(s_trial_nums)):
            t, n = np.unique(s_trial_nums, return_counts=True)
            repeats = t[n > 1]
            raise ValueError('Serial stream has repeated trial numbers: {0}.'.format(repeats))
        for tn in v_trial_nums:
            if tn not in s_trial_nums:
                logging.error('Trial number {0} found in Voyeur data file but not in serial stream!'.format(tn))

        trial_sets.append(v_trials)
        fieldsets.append(rec.Voyeur_data.Trials.colnames)
        sample_offsets.append(offset)
        offset += len(serial_st)
    sample_offsets.append(offset)  #append the last offset value to close the end of the session.
    grp._v_attrs['sample_rate_Hz'] = fs
    dest_file.root.events._v_attrs['sample_rate_Hz'] = fs
    # make superset of all fields. This will have fields that are not contained in all trials tb, so we'll have to
    # be careful for exemptions later on here.
    all_fields = []
    for fieldset in fieldsets:
        for field in fieldset:
            if not field in all_fields:
                all_fields.append(field)
    table_desc = dict()
    for i, f_name in enumerate(all_fields):
        for trs in trial_sets:
            if f_name in trs.colnames:
                try:
                    table_desc[f_name] = trs.coldescrs[f_name].copy(pos=i)
                except TypeError:  # if it's a string column, it needs some help.
                    # find the longest string that is here so that we can define the itemsize of the string.
                    max_len = 1  # this can't be 0.
                    for _ts in trial_sets:
                        for _t in _ts:
                            v = _t[f_name]
                            lv = len(v)
                            if lv > max_len:
                                max_len = lv
                    table_desc[f_name] = trs.coldescrs[f_name].copy(itemsize=max_len, pos=i)
            else:
                pass
    table_desc['run'] = tb.IntCol(pos=i+1)
    tbl = dest_file.create_table(where=grp,
                                 name=u'params_table',
                                 description=table_desc,
                                 filters=None,
                                 expectedrows=n_events)
    ev_time_array = []
    ev_off_time_array = []

    for rec in xrange(len(serial_dicts)):
        serial_dict = serial_dicts[rec]
        trials = trial_sets[rec]
        offset = sample_offsets[rec]
        trial_nums = trials[:]['trialNumber']
        tns = serial_dict.keys()
        tns.sort()
        for i in xrange(len(tns)):
            tn = tns[i]
            st_time = serial_dict[tn] + offset
            try:
                end_time = serial_dict[tns[i+1]] + offset
            except IndexError:  # last trial will cause this, catch it and use the end of the run from offsets.
                end_time = sample_offsets[rec+1]
            tr_id = np.where(trial_nums == tn)[0]
            if tr_id.size < 1:
                raise ValueError(u'Trial %i not found in voyeur trials table for run %i.'%(tn, rec))
            elif tr_id.size >1:
                raise ValueError(u'Multiple trials of number %i for run %i were found.' %(tn, rec))

            tr_id = tr_id[0]
            tr = trials[tr_id]
            row = tbl.row

            ev_time_array.append(st_time)
            ev_off_time_array.append(end_time)
            for field in all_fields:
                try:
                    row[field] = tr[field]
                except IndexError:
                    pass
            row['run'] = rec
            row.append()
    tbl.flush()
    evs = np.array([ev_time_array, ev_off_time_array]).T
    dest_file.create_carray(grp, u'events', obj=evs)
    dest_file.flush()
    return

def make_run_events(raw_kwd, dest_file):
    """

    :param raw_kwd:
    :param dest_file:
    :return:
    """
    assert isinstance(raw_kwd, tb.File)
    assert isinstance(dest_file, tb.File)
    grp = dest_file.create_group(u'/events',
                                 u'runs',
                                 title=u'start of voyeur run recordings',
                                 createparents=True)
    start = 0
    events = []
    end_events = []
    attr_sets = []
    for r in xrange(raw_kwd.root.recordings._v_nchildren):
        events.append(start)
        rec = raw_kwd.get_node(u'/recordings/{0:d}'.format(r))
        attr_sets.append(rec.Voyeur_data._v_attrs)
        start += len(rec.serial_trial)
        end_events.append(start)  # append the next start as the end of the current run.
    try:
        fs = rec.serial_trial.get_attr('sample_rate_Hz')
    except AttributeError:
        fs = 20833.
    grp._v_attrs['sample_rate_Hz'] = fs
    table_desc = {}
    for attrs in attr_sets:
        for k in attrs._v_attrnamesuser:
            if not table_desc.has_key(k):
                try:
                    a = tb.Atom.from_dtype(np.dtype((attrs[k].dtype, attrs[k].shape)))
                    table_desc[k] = tb.Col.from_atom(a)
                    # table_desc[k] = tb.Col.from_dtype((attrs[k].dtype, attrs[k].shape))
                except ValueError:  # string problems. string solutions
                    max_len = 1
                    for _a in attr_sets:
                        v = _a[k]
                        lv = len(v)
                        if lv > max_len:
                            max_len = lv
                    table_desc[k] = tb.StringCol(itemsize=max_len)
                except AttributeError:  # handle python objects (ie dict/tuple) by saving as a pickled string.
                    max_len = 1
                    for _a in attr_sets:
                        logging.debug('Saving Voyeur protocol attribute {0:s} as pickel object.'.format(k))
                        v = _a[k]
                        payload = pickle.dumps(v, protocol=0)  # dump as ascii (protocol 0).
                        lv = len(payload)
                        if lv > max_len:
                            max_len = lv
                    table_desc[k] = tb.StringCol(itemsize=max_len)
    tbl = dest_file.create_table(where=grp,
                                 name=u'params_table',
                                 description=table_desc,
                                 filters=None,
                                 expectedrows=raw_kwd.root.recordings._v_nchildren)
    for attrs in attr_sets:
        row = tbl.row
        for k in attrs._v_attrnamesuser:
            try:
                row[k] = attrs[k]
            except TypeError:  # handle python datatypes (ie dict/tuple) by saving as a pickled string.
                payload = pickle.dumps(attrs[k], protocol=0)  # dump as ascii (protocol 0).
                row[k] = payload
        row.append()
    tbl.flush()
    arr = np.array([events, end_events]).T
    dest_file.create_carray(grp, 'events', obj=arr)
    dest_file.flush()
    return


def get_trial_params(trial_events, trial_params, time):
    """
    Find the trial parameters that were uploaded immediately prior to a given time.

    :param trial_events: array of trial events (from /events/trials/events)
    :param trial_params: tb.Table for trial params_table (/events/trials/params_table)
    :param time: numeric time (in samples).
    :type trial_events: tb.Array
    :type trial_params: tb.Table
    :return:
    """

    tr_starts = trial_events[:, 0]
    tr_ends = trial_events[:, 1]
    try:
        trial_idx = np.where(tr_starts <= time)[0][-1]
        # print trial_idx
        tr = trial_params[trial_idx]
        assert tr_ends[trial_idx] > time, ('Strange, the time does not fall ',
                                           'within the end of the trial, something is wrong.')
    except IndexError:
        logging.error('No trial starts before time = %i, cannot retrieve trial information for this event.' % time)
        tr = None
    return tr


def _final_valve_handler(stream, stream_meta, trial):
    pass

def write_metadata(self, metadata_list, metadata_table):

    self.meta_table.append(metadata_list)
    self.meta_table.flush()

class GenericEventHandler(object):
    metadata_field_names = ()

    @property
    def stream_name(self):
        raise NotImplementedError('Event handler classes must define stream names to convert to events.')

    @property
    def metadata_field_names(self):
        raise NotImplementedError('Event handler class must define metadata_field_names or must specify None.')

    def __init__(self, raw_kwd, dest_file, *args, **kwargs):
        """
        Handles event creation from a trial. This means that it will convert events within a stream (ie TTL edges) as
        timestamps. With these timestamps, the class will bind metadata from a Voyeur H5 file to the events by matching
        the time when the event happened with the trial that is running at that time.

        To use for most purposes, inherit this class and overwrite the metadata_field_names tuple with the field names
        that are relevant to your stream.

        To record the voltage of the event (ie for an analog event like a laser pulse generator), add 'voltage' to the
        metadata_field_names tuple. (Find event times should be updated with a lower threshold, too.)

        :param raw_kwd:
        :param kwik:
        :param dest_file:
        :param args:
        :param kwargs:
        :return:
        """


        # TO BE OVERWRITTEN! This defines the field names (from the voyeur table)
        # that will be written as the event's metadata.

        assert isinstance(dest_file, tb.File)
        assert isinstance(raw_kwd, tb.File)
        self.dest_file = dest_file
        self.raw_kwd = raw_kwd

        stream = self.build_stream()
        self.stream = stream
        if stream is None:
            logging.warning('No "{0:s}" array found, not processing events.'.format(self.stream_name))
            return
        self._num_events = 0  # to use for estimation of table size. updated by find_event_times.
        self.events = self.find_event_times(stream)
        if self.events is None:
            return
        else:
            self.ev_grp = dest_file.create_group(u'/events', u'{0:s}'.format(self.stream_name), createparents=True)
            self.ev_grp._v_attrs['sample_rate_Hz'] = self.fs
            self.ev_events = dest_file.create_carray(self.ev_grp, 'events', obj=self.events)
            dest_file.flush()
            if self.metadata_field_names is not None:
                self.params_table = self.make_table()
                self.populate_params()
                dest_file.flush()
        return

    def build_stream(self):
        try:
            rec = 0
            nd_st = u'/recordings/{0:d}/{1:s}'.format(rec, self.stream_name)
            st_ex = self.raw_kwd.get_node(nd_st)
            st = np.array([], dtype=st_ex.dtype)
            for rec in xrange(self.raw_kwd.root.recordings._v_nchildren):
                nd_st = u'/recordings/{0:d}/{1:s}'.format(rec, self.stream_name)
                a = self.raw_kwd.get_node(nd_st)
                st = np.append(st, a.read())
            try:
                self.fs = a.get_attr('sample_rate_Hz')
            except AttributeError:
                self.fs = 20833.
        except tb.NoSuchNodeError:
                logging.warning(u'No {0:s} stream found for run {1:d} ({2:s})'.format(self.stream_name, rec, nd_st))
                st = None
        return st

    def find_event_times(self, stream):
        """

        :param stream:
        :return:
        """

        stream_log = stream > np.mean(stream)
        edges = np.convolve(stream_log, [1, -1])
        up_events = np.where(edges == 1)[0]
        down_events = np.where(edges == -1)[0]
        if up_events.size == 0:
            logging.warning(u'No events found in stream {0:s}'.format(self.stream_name))
            return None

        # assume that first up event should occur before the first down.
        if down_events[0] < up_events[0]:
            # if this isn't true, remove the first down event, as we probably don't have trial data for it.
            down_events = down_events[1:]
        # check if last down event was lost due to end of record. Append end of record time as last down time if so.
        if up_events.size - down_events.size == 1:
            down_events = np.append(down_events, stream.size)  # last down is after acquisition stops.
        assert up_events.size == down_events.size, u'Up and down events for stream {0:s} are not of equal size.'.format(
            self.stream_name)

        for i in xrange(1, len(down_events)):
            d = down_events[i]
            u = up_events[i]
            d2 = down_events[i-1]
            #if the previous down comes AFTER the last recorded up, we're missing an up event got problems.
            if (d - u) > (d - d2):
                raise ValueError(u'missing an up event in stream {0:s}.'.format(self.stream_name))
                # todo: try to recover from this.
        for i in xrange(len(up_events)-1):
            u = up_events[i]
            d = down_events[i]
            u2 = up_events[i+1]
            if (u2 - u) < (d - u):  # we're missing a down event if the next down comes after the next up.
                raise ValueError(u'missing a down event in stream {0:s}.'.format(self.stream_name))
                # todo: try to recover from this by adding fake down?

        events = np.array([up_events, down_events]).T
        self._num_events = events.shape[0]  #update estimate for making metadata table size.
        return events



    def make_table(self):

        # Make a table description dictionary for all of the fields in the metadata field names. This gets the column
        # description from the voyeur table, and replicates this column type for use in the event metadata table.
        #
        # Table column ordering (positions) will be consistent with the field ordering in metadata_field_names.

        #first get an example from the trial containing the first event:

        trial_params = self.dest_file.get_node(u'/events/trials/params_table')

        table_desc = dict()
        for i, f_name in enumerate(self.metadata_field_names):
            if f_name == 'voltage':  # we're going to record the voltage of the event (ie for laser), which we will measure from the stream, not take from the table, so we need to make a placeholder for it.
                table_desc['voltage'] = tb.IntCol(pos=i)  # int32 should be ok, but this can be flexiblized later..
            elif f_name in trial_params.colnames:
                try:
                    table_desc[f_name] = trial_params.coldescrs[f_name].copy(pos=i)
                except TypeError:  # if it's a string column, it needs some help.
                    # find the longest string that is here so that we can define the itemsize of the string.
                    max_len = 1
                    for _t in trial_params:
                        v = _t[f_name]
                        lv = len(v)
                        if lv > max_len:
                            max_len = lv
                    table_desc[f_name] = trial_params.coldescrs[f_name].copy(itemsize=max_len, pos=i)
            else:
                logging.warning('No field named "{0:s}" found in trials table!!!'.format(f_name))

        return self.dest_file.create_table(where=self.ev_grp, name='params_table', description=table_desc, filters=None,
                                           expectedrows=self._num_events)

    def populate_params(self):
        """
        This populates the event parameters with parameters corresponding to the trial that is active at the time of the
        event onset.

        THIS can be overwritten for special cases (ie where multiple events within a trial have different parameters,
        but it will take some effort.

        :return:
        """
        tr_events = self.dest_file.get_node('/events/trials/events')
        tr_params = self.dest_file.get_node('/events/trials/params_table')
        table = self.params_table
        for event in self.events:
            row = self.params_table.row
            tr = get_trial_params(tr_events, tr_params, event[0])
            if tr:
                for field in self.metadata_field_names:
                    if field.lower() == 'voltage':
                        row[field] = self.read_stream_voltage(event[0])
                    row[field] = tr[field]
            row.append()
        table.flush()


    def handle_metadata(self):

        # make a list of tuples of metadata.
        all_md = []
        for event in self.on_events:
            # find the proximal early trial start. This will have the metadata regarding this event.
            tr_start_samples = np.array(self.trial_starts.values())
            tr_start_samples.sort()
            tr_start_sample = tr_start_samples[tr_start_samples < self.event[0]][-1]
            ev_tr_num = self.trial_starts[tr_start_sample]
            tr_nums = self.table['trialNumber']
            tr_idx = np.where(tr_nums == ev_tr_num)[0]
            if len(tr_idx) > 1:
                raise ValueError(u'WARNING multiple trials with same numbering. This will result in ambiguous metadata.')
            tr_idx = tr_idx[0]
            ev_tr = self.table[tr_idx]
            md = []
            for field in self.metadata_field_names:  # ordering of metadata_field_names tuple/list is constant so we can use it for the ordering of the list to append.
                if field == 'voltage':
                    md.append(self.read_stream_voltage(event))
                else:
                    md.append(ev_tr[field])
            all_md.append(tuple(md))

        return


    def read_stream_voltage(self, event_on_time):
        """

        :param event_on_time: sample where to read the stream to get the voltage reading.
        :return:
        """

        return self.stream[event_on_time]



class FinalValveEventHandler(GenericEventHandler):
    """
    Handles odor information. Since every odor presentation involves a final valve opening, this is a great event in which
    to store odor info. This should be essentially a default event handler, all we have to do is define metadata field
    names.

    The super init will do everything else for us.
    """

    metadata_field_names = ('odor',
                            'odorconc',
                            'NitrogenFlow_1',
                            'dillution',
                            'vialconc',
                            'fvdur',
                            'fvOnTime')
    stream_name = u'finalvalve'

class FinalValveEventHandlerArduino(GenericEventHandler):
    """
    This handles final valve time handling for files that do not have the FV trig recorded (mouse 3861) by using the
    time of the Arduino recorded trial parameters recieval and the Arduino reported FV time.

    This will not impact recordings which have the FV recorded within their raw.kwd file, as it will detect that these
    are present and return without doing anything.
    """

    metadata_field_names = ('odor', 'odorconc', 'NitrogenFlow_1', 'dillution', 'vialconc', 'fvdur', 'fvOnTime')
    stream_name = u'finalvalve'

    def __init__(self, raw_kwd, dest_file, *args, **kwargs):
        assert isinstance(dest_file, tb.File)
        if self.check_no_stream(raw_kwd):
            logging.info('No stream found for {0:s}, creating events from Arduino record.'.format(self.stream_name))
            super(FinalValveEventHandlerArduino, self).__init__(raw_kwd, dest_file, *args, **kwargs)
            self.ev_grp._f_setattr('Warning', 'Events processed using FinalValveEventHandlerArduino!')
        else:
            logging.info('Stream found for {0:s}, NOT creating events from Arduino record.')
            return

    def check_no_stream(self, raw_kwd):
        """
        Returns true of stream of self.stream_name is not present in the raw.kwd.

        :param raw_kwd:
        :return:
        """
        try:
            raw_kwd.get_node('/recordings/0/{0:s}'.format(self.stream_name))
            return False
        except tb.NoSuchNodeError:
            return True

    def build_stream(self):
        self.fs = 20833.
        return []

    def find_event_times(self, stream):
        trial_events = self.dest_file.get_node(u'/events/trials/events')
        trial_params = self.dest_file.get_node(u'/events/trials/params_table')

        try:
            recording = self.raw_kwd.get_node('/recordings/0/neural')
            fs = recording._v_attrs['sample_rate_Hz']
        except:
            logging.warning('No sampling rate found in node raw.kwd/recordings/0/neural, '
                            'defaulting to 20833. Hz (whisper).')
            fs = 20833.

        up_events = []
        down_events = []
        for i in xrange(len(trial_events)):
            tev = trial_events[i]
            tparam = trial_params[i]

            if tparam['fvOnTime']:  # zero here means no FV opened.
                t_start_acq = tev[0] # when serial was recieved by acquistion system
                t_start_arduino = tparam['paramsgottime']  # just before serial was printed in arduino time
                t_fv_on_arduino = tparam['fvOnTime']
                fv_on_diff_ms = t_fv_on_arduino - t_start_arduino  # diff in time between trial start and fv on IN MS.
                fv_on_diff_samp = np.int(fs/1000. * fv_on_diff_ms)
                fv_on_acq = t_start_acq + fv_on_diff_samp
                fv_off_ms = tparam['fvdur'] + fv_on_diff_ms
                fv_off_samp = np.int(fs/1000. * fv_off_ms)
                fv_off_acq = fv_on_acq + fv_off_samp
                up_events.append(fv_on_acq)
                down_events.append(fv_off_acq)
        down_events = np.array(down_events, dtype=np.int)
        up_events = np.array(up_events, dtype=np.int)
        events = np.array([up_events, down_events]).T
        self._num_events = events.shape[0]
        return events


class LickEventHandler(GenericEventHandler):
    def __init__(self, *args, **kwargs):
        pass


class LaserPulseEventHandler(GenericEventHandler):
    def __init__(self, *args, **kwargs):
        pass

class WaterEventHandler(GenericEventHandler):
    def __init__(self, *args, **kwargs):
        pass


handlers = [FinalValveEventHandler, LaserPulseEventHandler, WaterEventHandler, FinalValveEventHandlerArduino]