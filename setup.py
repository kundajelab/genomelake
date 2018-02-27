from setuptools.extension import Extension
from setuptools import setup, find_packages
from codecs import open
from os import path

from Cython.Build import cythonize

setup(
    name='genomelake',

    version='0.1.1',

    description='Simple and efficient random access to genomic data for deep learning models.',
    long_description='',

    url='https://github.com/kundajelab/genomelake',

    author='Chuan-Sheng Foo',
    author_email='csfoo@cs.stanford.edu',

    license='BSD-3',

    classifiers=[
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Science/Research',
        'Intended Audience :: Developers',
        'Topic :: Scientific/Engineering',
        'License :: OSI Approved :: BSD License',
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: 3.6',
        'Operating System :: POSIX :: Linux',
    ],

    keywords='deeplearning neuralnets genomics',

    packages=find_packages(exclude=['contrib', 'docs', 'tests*']),

    ext_modules=cythonize([Extension('genomelake.util',
                                     ['genomelake/util.pyx'])]),

    install_requires=['bcolz>=1.1', 'numpy', 'pybedtools',
                      'pysam'],
)
