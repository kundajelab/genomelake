from __future__ import absolute_import, division, print_function
import json
import numpy as np
import os
import six

import bcolz
from pybedtools import BedTool
import pyBigWig
from pysam import FastaFile

from .util import makedirs
from .util import one_hot_encode_sequence
from .util import nan_to_zero

NUM_SEQ_CHARS = 4

_blosc_params = bcolz.cparams(clevel=5, shuffle=bcolz.SHUFFLE, cname='lz4')

_array_writer = {
    'numpy': lambda arr, path: np.save(path, arr),
    'bcolz': lambda arr, path: bcolz.carray(
        arr, rootdir=path, cparams=_blosc_params, mode='w').flush()
}


def extract_fasta_to_file(fasta, output_dir, mode='bcolz', overwrite=False):
    assert mode in _array_writer

    makedirs(output_dir, exist_ok=overwrite)
    fasta_file = FastaFile(fasta)
    file_shapes = {}
    for chrom, size in zip(fasta_file.references, fasta_file.lengths):
        data = np.empty((size, NUM_SEQ_CHARS), dtype=np.float32)
        seq = fasta_file.fetch(chrom)
        one_hot_encode_sequence(seq, data)
        file_shapes[chrom] = data.shape
        _array_writer[mode](data, os.path.join(output_dir, chrom))

    with open(os.path.join(output_dir, 'metadata.json'), 'w') as fp:
        json.dump({'file_shapes': file_shapes,
                   'type': 'array_{}'.format(mode),
                   'source': fasta}, fp)


def extract_bigwig_to_file(bigwig, output_dir, mode='bcolz', dtype=np.float32,
                           overwrite=False, nan_as_zero=True):
    assert mode in _array_writer

    makedirs(output_dir, exist_ok=overwrite)
    bw = pyBigWig.open(bigwig)
    chrom_sizes = bw.chroms()
    file_shapes = {}
    for chrom, size in six.iteritems(chrom_sizes):
        data = np.empty(size)
        data = bw.values(chrom, 0, size, numpy=True)
        if nan_as_zero:
            nan_to_zero(data)
        _array_writer[mode](data.astype(dtype),
                            os.path.join(output_dir, chrom))
        file_shapes[chrom] = data.shape
    bw.close()

    with open(os.path.join(output_dir, 'metadata.json'), 'w') as fp:
        json.dump({'file_shapes': file_shapes,
                   'type': 'array_{}'.format(mode),
                   'source': bigwig}, fp)


def read_genome_sizes(genome_file):
    with open(genome_file) as fp:
        chr2size = {}
        for line in fp:
            chrom, size = line.split()
            chr2size[chrom] = int(size)
    return chr2size


def load_directory(base_dir, in_memory=False):
    with open(os.path.join(base_dir, 'metadata.json'), 'r') as fp:
        metadata = json.load(fp)

    if metadata['type'] == 'array_numpy':
        mmap_mode = None if in_memory else 'r'
        data = {chrom: np.load('{}.npy'.format(os.path.join(base_dir, chrom)),
                                mmap_mode=mmap_mode)
                for chrom in metadata['file_shapes']}

        for chrom, shape in six.iteritems(metadata['file_shapes']):
            if data[chrom].shape != tuple(shape):
                raise ValueError('Inconsistent shape found in metadata file: '
                                 '{} - {} vs {}'.format(chrom, shape,
                                                        data[chrom].shape))
    elif metadata['type'] == 'array_bcolz':
        data = {chrom: bcolz.open(os.path.join(base_dir, chrom), mode='r')
                for chrom in metadata['file_shapes']}
        if in_memory:
            data = {k: data[k].copy() for k in data.keys()}

        for chrom, shape in six.iteritems(metadata['file_shapes']):
            if data[chrom].shape != tuple(shape):
                raise ValueError('Inconsistent shape found in metadata file: '
                                 '{} - {} vs {}'.format(chrom, shape,
                                                        data[chrom].shape))
    else:
        raise ValueError('Can only extract from array_bcolz and array_numpy')

    return data
