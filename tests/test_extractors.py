from genomelake.extractors import FastaExtractor
import numpy as np
from pybedtools import Interval
import pytest

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
