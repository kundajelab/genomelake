# genomelake
[![CircleCI](https://circleci.com/gh/kundajelab/genomelake.svg?style=svg)](https://circleci.com/gh/kundajelab/genomelake)[![Coverage Status](https://coveralls.io/repos/github/kundajelab/genomelake/badge.svg)](https://coveralls.io/github/kundajelab/genomelake)

Efficient random access to genomic data for deep learning models.

Supports the following types of input data:

- bigwig
- DNA sequence

genomelake extracts signal from genomic inputs in provided BED intervals.

## Requirements
- python 2.7 or 3.5
- bcolz
- cython
- numpy
- pybedtools
- pysam

## Installation
Clone the repository and run:

`python setup.py install`

## Getting started: training a protein-DNA binding model
Extract genome-wide sequence data into a genomelake data source:
```python
from genomelake.backend import extract_fasta_to_file

genome_fasta = "/mnt/data/annotations/by_release/hg19.GRCh37/hg19.genome.fa"
genome_data_directory = "./hg19_data_directory"
extract_fasta_to_file(genome_fasta, genome_data_directory)
```

Using a BED intervals file with labels, a genome data source, and genomelake's `ArrayExtractor`, generate input DNA sequences and labels:
```python
import pybedtools
from genomelake.extractors import ArrayExtractor
import numpy as np

def batch_iter(iterable, batch_size):
    it = iter(iterable)
    try:
        while True:
            values = []
            for n in range(batch_size):
                values += (next(it),)
            yield values
    except StopIteration:
        yield values

def generate_inputs_and_labels(intervals_file, data_source, batch_size=128):
    bt = pybedtools.BedTool(intervals_file)
    extractor = ArrayExtractor(data_source)
    intervals_generator = batch_iter(bt, batch_size)
    for intervals_batch in intervals_generator:
    	inputs = extractor(intervals_batch)
	labels = []
	for interval in intervals_batch:
	    labels.append(float(interval.name))
        labels = np.array(labels)
        yield inputs, labels
```

Train a keras model of JUND binding to DNA using 101 base pair intervals and labels in `./examples/JUND.HepG2.chr22.101bp_intervals.tsv.gz`:
```python
from keras.models import Sequential
from keras.layers import Conv1D, Flatten, Dense

intervals_file = "./examples/JUND.HepG2.chr22.101bp_intervals.tsv.gz"
inputs_labels_generator = generate_inputs_and_labels(intervals_file, genome_data_directory)

model = Sequential()
model.add(Conv1D(15, 25, input_shape=(101, 4)))
model.add(Flatten())
model.add(Dense(1, activation='sigmoid'))

model.compile(loss='binary_crossentropy', optimizer='adam', metrics=['accuracy'])
model.fit_generator(inputs_labels_generator, steps_per_epoch=100)
```

Here is the expected result:
```
100/100 [==============================] - 7s - loss: 0.0584 - acc: 0.9905 
```

## License
genomelake is released under the BSD-3 license. See ``LICENSE`` for details.
