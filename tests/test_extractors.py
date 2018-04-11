from genomelake import backend
from genomelake.extractors import ArrayExtractor, FastaExtractor
import numpy as np
from pybedtools import Interval
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
