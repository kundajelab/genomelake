import os
from setuptools.extension import Extension
from setuptools import setup, find_packages

def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()

setup(
    name='genomelake',

    version='0.1.3',

    description='Simple and efficient random access to genomic data for deep learning models.',
    long_description=read('README.md'),

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

    setup_requires=['cython'],

    ext_modules=[Extension('genomelake.util', ['genomelake/util.pyx'])],

    install_requires=['bcolz>=1.1', 'numpy', 'pybedtools',
                      'pyBigWig>=0.3.2', 'pysam', 'six>=1.9.0'],
)
