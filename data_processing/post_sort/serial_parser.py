__author__ = 'chris'

import numpy as np


def parse_serial_stream(serial_stream, fs=20833., word_len=2, baudrate=300):
    """

    :param serial_stream: array containing serial stream.
    :param fs: sampling frequency of serial stream
    :param word_len: number of bytes per word (ie 16 bit int is 2 bytes)
    :param baudrate: baudrate of serial transmission.
    :return:
    """

    if serial_stream.ndim == 2:
        try:
            shortdim = serial_stream.shape.find(1)
            if shortdim == 1:
                serial_stream = serial_stream[:, 0]
            else:
                serial_stream = serial_stream[0, :]
        except ValueError:
            raise ValueError('parse_serial_stream input must be a 1d array or 2d with dimension of length 1.')

    start_bit = 1
    end_bit = 1
    bit_len = float(fs)/baudrate

    log_stream = (serial_stream > np.max(serial_stream)/2)
    edges = np.convolve(log_stream, [1,-1])
    ups = np.where(edges == 1)[0]
    downs = np.where(edges == -1)[0]

    #check to see if baudrate/fs combo is reasonable by looking for bits of the length specified by them in the stream.
    diff = (downs-ups) / bit_len
    cl = 0
    for i in range(1,4):  # look for bit widths of 1,2,3 consecutive.
        diff -= 1.
        g = diff > -.1
        l = diff < .1
        cl += np.sum(g*l)
    if cl < 4:
        print 'WARNING: bits do not appear at timings consistent with selected \nsampling frequency and baud rate ' \
              'combination: (fs: %i, baud: %i). This warning may be erroneous for short recordings.' % (fs, baudrate)

    start_id = downs-ups > (bit_len * word_len * (8 + start_bit + end_bit))
    start_samples = downs[start_id]

    word_bits = bit_len * word_len * (8+start_bit+end_bit)
    bit_sample_spacing = np.linspace(bit_len/2, word_bits - bit_len/2, word_len * 10)
    bit_sample_spacing = bit_sample_spacing.round().astype(int)

    _bits_tmp = np.zeros(bit_sample_spacing.shape, dtype=np.bool)
    _bytes = np.empty((word_len, 8), dtype=bool)

    trial_time_dict = {}
    for start in start_samples[:-1]:
    #     print start
        bit_samples = bit_sample_spacing + start
        bits = log_stream[bit_samples]
        for i in range(word_len):
            _bytes[i,:] = bits[i*9+1:(i+1)*9]
        number = _bytes_to_int(_bytes)
        trial_time_dict[number] = start
    return trial_time_dict


def _bytes_to_int(bytes, endianness='little'):
    """

    :param bytes: np boolean array with shape == (number_bytes, 8)
    :param endianness of the encoded bytes
    :return:
    """

    bits = bytes.ravel()
    if endianness == 'big':
        bits = np.flipud(bits)
    num = 0
    for e, bit in enumerate(bits):
        num += bit * 2 ** e
    return num