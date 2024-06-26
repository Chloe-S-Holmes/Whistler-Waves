import struct
import collections
import h5py
import numpy as np

class DataSet:
    WAVEDESC = collections.namedtuple('WAVEDESC',
                ['descriptor_name', 'template_name', 'comm_type', 'comm_order',
                 'wave_descriptor', 'user_text', 'res_desc1', 'trigtime_array', 'ris_time_array',
                 'res_array1', 'wave_array_1', 'wave_array_2', 'res_array2', 'res_array3',
                 'instrument_name', 'instrument_number', 'trace_label', 'reserved1', 'reserved2',
                 'wave_array_count', 'pnts_per_screen', 'first_valid_pnt', 'last_valid_pnt',
                 'first_point', 'sparsing_factor', 'segment_index', 'subarray_count', 'sweeps_per_acq',
                 'points_per_pair', 'pair_offset', 'vertical_gain', 'vertical_offset', 'max_value',
                 'min_value', 'nominal_bits', 'nom_subarray_count', 'horiz_interval', 'horiz_offset',
                 'pixel_offset', 'vertunit', 'horunit', 'horiz_uncertainty',
                 'tt_second', 'tt_minute', 'tt_hours', 'tt_days', 'tt_months', 'tt_year', 'tt_unused',
                 'acq_duration', 'record_type', 'processing_done', 'reserved5', 'ris_sweeps',
                 'timebase', 'vert_coupling', 'probe_att', 'fixed_vert_gain', 'bandwidth_limit',
                 'vertical_vernier', 'acq_vert_offset', 'wave_source'])
    
    RECORD_TYPES = ['single_sweep', 'interleaved', 'histogram', 'graph', 'filter_coefficient',
                    'complex', 'extrema', 'sequence_obsolete', 'centered_RIS', 'peak_detect']
    
    PROCESSING_TYPES = ['no_processing', 'fir_filter', 'interpolated', 'sparsed',
                        'autoscaled', 'no_result', 'rolling', 'cumulative']
    
    TIMEBASE_IDS = ['1 ps', '2 ps', '5 ps', '10 ps', '20 ps', '50 ps', '100 ps', '200 ps', '500 ps',
                    '1 ns', '2 ns', '5 ns', '10 ns', '20 ns', '50 ns', '100 ns', '200 ns', '500 ns',
                    '1 us', '2 us', '5 us', '10 us', '20 us', '50 us', '100 us', '200 us', '500 us',
                    '1 ms', '2 ms', '5 ms', '10 ms', '20 ms', '50 ms', '100 ms', '200 ms', '500 ms',
                    '1 s',  '2 s',  '5 s',  '10 s',  '20 s',  '50 s',  '100 s',  '200 s',  '500 s',
                    '1 ks', '2 ks', '5 ks']   # these are per division; ALSO: 100 corresponds to EXTERNAL
    
    VERT_GAIN_IDS = ['1 uV', '2 uV', '5 uV', '10 uV', '20 uV', '50 uV', '100 uV', '200 uV', '500 uV',
                     '1 mV', '2 mV', '5 mV', '10 mV', '20 mV', '50 mV', '100 mV', '200 mV', '500 mV',
                     '1 V',  '2 V',  '5 V',  '10 V',  '20 V',  '50 V',  '100 V',  '200 V',  '500 V',
                     '1 kV', '2 kV', '5 kV', '10 kV']   # these are per division; pp added the last 3
    
    VERT_COUPLINGS = ['DC 50 Ohms', 'ground', 'DC 1 MOhm', 'ground', 'AC 1 MOhm']
    
    def __init__(self, fname):
        self.fname = fname

    def print_tree(self):
        def printname(name):
            print(name)
        with h5py.File(self.fname, 'r') as f:
             f.visit(printname)
    
    def read_data(self, chan):
        with h5py.File(self.fname, 'r') as f:
            c_data = np.asarray(f[f'Acquisition/LeCroy_scope/Channel{int(chan)}'])
        return c_data

    def read_header(self, chan, silent = False):
        with h5py.File(self.fname, 'r') as f:
            c_header = np.asarray(f[f'Acquisition/LeCroy_scope/Headers/Channel{int(chan)}'])
        if not silent:
            self.check_header(c_header)

    def read_time(self):
        with h5py.File(self.fname, 'r') as f:
            t_array = np.asarray(f[f'Acquisition/LeCroy_scope/time'])
        return t_array

    def read_positions(self):
        with h5py.File(self.fname, 'r') as f:
            t_array = np.asarray(f[f'Control/Positions/positions_setup_array'])
        return t_array
    
    def check_header(self, a):
        hdr = []
        for i, val in enumerate(struct.unpack('=16s16shhllllllllll16sl16shhlllllllllhhffffhhfdd48s48sfdBBBBhhfhhhhhhfhhffh', a[0])):
            if isinstance(val, bytes):
                val = val.decode('ascii').strip("\x00")
            # if i < len(_trc_description_fields):
                #print(_trc_description_fields[i],val)
            hdr.append(val)
        header = self.WAVEDESC(*hdr)
        print(header)
        print('Record Type: ', self.RECORD_TYPES[header.record_type], '\n',
          'Processing Type: ', self.PROCESSING_TYPES[header.processing_done], '\n',
          'Time/div: ', self.TIMEBASE_IDS[header.timebase], '\n',
          'V/div: ', self.VERT_GAIN_IDS[header.fixed_vert_gain], '\n',
          'Input Coupling: ', self.VERT_COUPLINGS[header.vert_coupling])