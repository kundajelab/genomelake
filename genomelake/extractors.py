from __future__ import absolute_import, division, print_function
import numpy as np
import six

import bcolz
from pybedtools import BedTool
from pybedtools import Interval
from pysam import FastaFile

from . import backend
from .util import one_hot_encode_sequence

NUM_SEQ_CHARS = 4


class BaseExtractor(object):
    dtype = np.float32

    def __init__(self, datafile, **kwargs):
        self._datafile = datafile

    def __call__(self, intervals, out=None, **kwargs):
        data = self._check_or_create_output_array(intervals, out)
        self._extract(intervals, data, **kwargs)
        return data

    def _check_or_create_output_array(self, intervals, out):
        width = intervals[0].stop - intervals[0].start
        output_shape = self._get_output_shape(len(intervals), width)

        if out is None:
            out = np.zeros(output_shape, dtype=self.dtype)
        else:
            if out.shape != output_shape:
                raise ValueError('out array has incorrect shape: {} '
                                 '(need {})'.format(out.shape, output_shape))
            if out.dtype != self.dtype:
                raise ValueError('out array has incorrect dtype: {} '
                                 '(need {})'.format(out.dtype, self.dtype))
        return out

    def _extract(self, intervals, out, **kwargs):
        'Subclassses should implement this and return the data'
        raise NotImplementedError

    @staticmethod
    def _get_output_shape(num_intervals, width):
        'Subclasses should implement this and return the shape of output'
        raise NotImplementedError


class ArrayExtractor(BaseExtractor):

    def __init__(self, datafile, in_memory=False, **kwargs):
        super(ArrayExtractor, self).__init__(datafile, **kwargs)
        self._data = backend.load_directory(datafile, in_memory=in_memory)
        self.multiprocessing_safe = in_memory

        arr = next(iter(self._data.values()))
        # The reason why we do this is because bcolz doesn't support ellipsis
        # in indexing, unlike numpy.
        if isinstance(arr, bcolz.carray):
            def _mm_extract(self, intervals, out, **kwargs):
                mm_data = self._data
                for index, interval in enumerate(intervals):
                    out[index] = mm_data[interval.chrom][interval.start:interval.stop]
        else:
            def _mm_extract(self, intervals, out, **kwargs):
                mm_data = self._data
                for index, interval in enumerate(intervals):
                    out[index] = mm_data[interval.chrom][..., interval.start:interval.stop]

        # output shape method
        shape = arr.shape
        if len(shape) == 1:
            def _get_output_shape(num_intervals, width):
                return (num_intervals, width)
        elif len(shape) == 2:
            def _get_output_shape(num_intervals, width):
                return (num_intervals, width, shape[1])
        else:
            raise ValueError('Can only extract from 1D/2D arrays')

        self._mm_extract = _mm_extract.__get__(self)
        self._extract = self._mm_extract
        self._get_output_shape = staticmethod(_get_output_shape).__get__(self)


class FastaExtractor(BaseExtractor):

    def __init__(self, datafile, use_strand=False, **kwargs):
        super(FastaExtractor, self).__init__(datafile, **kwargs)
        self.use_strand = use_strand

    def _extract(self, intervals, out, **kwargs):
        fasta = FastaFile(self._datafile)

        for index, interval in enumerate(intervals):
            seq = fasta.fetch(str(interval.chrom), interval.start,
                              interval.stop)
            one_hot_encode_sequence(seq, out[index, :, :])

            # reverse-complement seq the negative strand
            if self.use_strand and interval.strand == "-":
                out[index, :, :] = out[index, ::-1, ::-1]

        return out

    @staticmethod
    def _get_output_shape(num_intervals, width):
        return (num_intervals, width, NUM_SEQ_CHARS)


class BigwigExtractor(BaseExtractor):

    def __init__(self, datafile, **kwargs):
        super(BigwigExtractor, self).__init__(datafile, **kwargs)
        self._verbose = kwargs.get('verbose', False)

    def _extract(self, intervals, out, **kwargs):
        out[:] = self._bigwig_extractor(self._datafile, intervals,
                                        **kwargs)

        return out

    @staticmethod
    def _get_output_shape(num_intervals, width):
        return (num_intervals, width)

    @staticmethod
    def _bigwig_extractor(datafile, intervals, out=None, **kwargs):
        if out is None:
            width = intervals[0].stop - intervals[0].start
            out = np.zeros((len(intervals), width))

        bw = pyBigWig.open(datafile)
        for index, interval in enumerate(intervals):
            out[index] = bw.values(interval.chrom, interval.start,
                                   interval.stop, numpy=True)
        bw.close()

        return out
