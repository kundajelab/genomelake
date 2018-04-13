from genomelake import backend
from genomelake.extractors import ArrayExtractor, BigwigExtractor, FastaExtractor
import numpy as np
from pybedtools import Interval
import pyBigWig
import pytest

array_extractor_fasta_params = [("numpy", True),
                                ("numpy", False),
                                ("bcolz", True),
                                ("bcolz", False)]

def test_fasta_extractor_valid_intervals():
    extractor = FastaExtractor('tests/data/fasta_test.fa')
    intervals = [Interval('chr1', 0, 10),
                 Interval('chr2', 0, 10)]
    expected_data = np.array(
        [[[ 1.  ,  0.  ,  0.  ,  0.  ],
          [ 0.  ,  1.  ,  0.  ,  0.  ],
          [ 0.  ,  1.  ,  0.  ,  0.  ],
          [ 0.  ,  0.  ,  1.  ,  0.  ],
          [ 0.  ,  0.  ,  0.  ,  1.  ],
          [ 1.  ,  0.  ,  0.  ,  0.  ],
          [ 0.  ,  1.  ,  0.  ,  0.  ],
          [ 0.  ,  1.  ,  0.  ,  0.  ],
          [ 0.  ,  0.  ,  1.  ,  0.  ],
          [ 0.  ,  0.  ,  0.  ,  1.  ]],

         [[ 1.  ,  0.  ,  0.  ,  0.  ],
          [ 0.  ,  1.  ,  0.  ,  0.  ],
          [ 0.  ,  0.  ,  1.  ,  0.  ],
          [ 0.  ,  0.  ,  0.  ,  1.  ],
          [ 0.25,  0.25,  0.25,  0.25],
          [ 1.  ,  0.  ,  0.  ,  0.  ],
          [ 0.  ,  1.  ,  0.  ,  0.  ],
          [ 0.  ,  0.  ,  1.  ,  0.  ],
          [ 0.  ,  0.  ,  0.  ,  1.  ],
          [ 0.25,  0.25,  0.25,  0.25]]], dtype=np.float32)
    data = extractor(intervals)
    assert (data == expected_data).all()


def test_fasta_extractor_over_chr_end():
    extractor = FastaExtractor('tests/data/fasta_test.fa')
    intervals = [Interval('chr1', 0, 100),
                 Interval('chr1', 1, 101)]
    with pytest.raises(ValueError):
        data = extractor(intervals)

@pytest.mark.parametrize("mode,in_memory", array_extractor_fasta_params)
def test_array_extractor_fasta(mode, in_memory):
    data_dir = 'tests/data/fasta_test_dir_{}_{}'.format(mode, in_memory)
    backend.extract_fasta_to_file(
        'tests/data/fasta_test.fa',
        data_dir,
        mode=mode,
        overwrite=True)
    extractor = ArrayExtractor(data_dir, in_memory=in_memory)
    intervals = [Interval('chr1', 0, 10),
                 Interval('chr2', 0, 10)]
    expected_data = np.array(
        [[[ 1.  ,  0.  ,  0.  ,  0.  ],
          [ 0.  ,  1.  ,  0.  ,  0.  ],
          [ 0.  ,  1.  ,  0.  ,  0.  ],
          [ 0.  ,  0.  ,  1.  ,  0.  ],
          [ 0.  ,  0.  ,  0.  ,  1.  ],
          [ 1.  ,  0.  ,  0.  ,  0.  ],
          [ 0.  ,  1.  ,  0.  ,  0.  ],
          [ 0.  ,  1.  ,  0.  ,  0.  ],
          [ 0.  ,  0.  ,  1.  ,  0.  ],
          [ 0.  ,  0.  ,  0.  ,  1.  ]],

         [[ 1.  ,  0.  ,  0.  ,  0.  ],
          [ 0.  ,  1.  ,  0.  ,  0.  ],
          [ 0.  ,  0.  ,  1.  ,  0.  ],
          [ 0.  ,  0.  ,  0.  ,  1.  ],
          [ 0.25,  0.25,  0.25,  0.25],
          [ 1.  ,  0.  ,  0.  ,  0.  ],
          [ 0.  ,  1.  ,  0.  ,  0.  ],
          [ 0.  ,  0.  ,  1.  ,  0.  ],
          [ 0.  ,  0.  ,  0.  ,  1.  ],
          [ 0.25,  0.25,  0.25,  0.25]]], dtype=np.float32)
    data = extractor(intervals)
    assert (data == expected_data).all()

@pytest.fixture
def test_bigwig_and_intervals():
    bw_path = "tests/data/test_bigwig.bw"
    intervals = [Interval('chr1', 0, 10),
                 Interval('chr2', 0, 10)]
    expected_chr1 = np.array([0.1] * 10, dtype=np.float32)
    expected_chr2 = np.array([0] + [9]*9, dtype=np.float32)
    expected_data = np.stack([expected_chr1, expected_chr2])

    return (bw_path, intervals, expected_data)

@pytest.mark.parametrize("mode,in_memory", array_extractor_fasta_params)
def test_array_extractor_bigwig(test_bigwig_and_intervals, mode, in_memory):
    bw_path, intervals, expected_data = test_bigwig_and_intervals
    bw_dir_path = "{}.dir".format(bw_path)
    backend.extract_bigwig_to_file(
        bw_path, bw_dir_path, mode=mode, overwrite=True)
    extractor = ArrayExtractor(bw_dir_path, in_memory=in_memory)

    data = extractor(intervals)
    assert (data == expected_data).all()


def test_bigwig_extractor(test_bigwig_and_intervals):
    bw_path, intervals, expected_data = test_bigwig_and_intervals
    extractor = BigwigExtractor(bw_path)
    data = extractor(intervals)
    assert (data == expected_data).all()
