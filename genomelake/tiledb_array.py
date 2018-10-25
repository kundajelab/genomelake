from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

import os.path
import shutil

import numpy as np
import tiledb

GENOME_DOMAIN_NAME = "genome_coord"
SECONDARY_DOMAIN_NAME = "signal_coord"
GENOME_VALUE_NAME = "v"

DEFAULT_GENOME_TILE_EXTENT = 9000
DEFAULT_COMPRESSOR = "blosc-lz"
DEFAULT_COMPRESSOR_LEVEL = -1


def write_tiledb(arr, path, overwrite=True):
    """Write a tiledb to disk.
    """
    if os.path.exists(path) and os.path.isdir(path) and overwrite:
        shutil.rmtree(path)

    if os.path.exists(path):
        raise FileExistsError("Output path {} already exists".format(path))

    ctx = tiledb.Ctx()

    n = arr.shape[0]
    n_tile_extent = min(DEFAULT_GENOME_TILE_EXTENT, n)

    d1 = tiledb.Dim(
        ctx, GENOME_DOMAIN_NAME, domain=(0, n - 1), tile=n_tile_extent, dtype="uint32"
    )

    if arr.ndim == 1:
        domain = tiledb.Domain(ctx, d1)

    elif arr.ndim == 2:
        m = arr.shape[1]
        d2 = tiledb.Dim(
            ctx, SECONDARY_DOMAIN_NAME, domain=(0, m - 1), tile=m, dtype="uint32"
        )
        domain = tiledb.Domain(ctx, d1, d2)

    else:
        raise ValueError("tiledb backend only supports 1D or 2D arrays")

    v = tiledb.Attr(
        ctx,
        GENOME_VALUE_NAME,
        compressor=(DEFAULT_COMPRESSOR, DEFAULT_COMPRESSOR_LEVEL),
        dtype="float32",
    )

    schema = tiledb.ArraySchema(
        ctx, domain=domain, attrs=(v,), cell_order="row-major", tile_order="row-major"
    )
    A = tiledb.DenseArray.create(path, schema)

    values = arr.astype(np.float32)

    with tiledb.DenseArray(ctx, path, mode="w") as A:
        A[:] = {GENOME_VALUE_NAME: values}


def load_tiledb(path):
    return TDBDenseArray(path)


class TDBDenseArray(object):
    """A read-only wrapper of tiledb.DenseArray"""

    def __init__(self, array):
        if isinstance(array, tiledb.DenseArray):
            self._arr = array
        else:
            ctx = tiledb.Ctx()
            self._arr = tiledb.DenseArray(ctx, array, mode="r")

    def __getitem__(self, key):
        return self._arr[key][GENOME_VALUE_NAME]

    def __setitem__(self, key, item):
        raise NotImplemented("TDBDenseArray is read-only")

    @property
    def shape(self):
        return self._arr.shape

    @property
    def ndim(self):
        return self._arr.ndim
