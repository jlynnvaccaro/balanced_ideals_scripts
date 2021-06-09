"""Open file with automatic (de)compression based on filename"""
# David Dumas <david@dumas.io>
# This work is placed in the public domain.
import os


def zopen(fn, *args, **kwargs):
    """Open a file, transparently decompressing on read or compressing on write
    if the file extension indicates bzip2 compression (.bz2) or gzip compression
    (.gz)."""
    # TODO: Handle (kw)args only supported by one of the backends gracefully
    # TODO: gzip, bzip use binary mode by default, open() uses text.  Should
    # something be done to enforce a consistent default?
    _, ext = os.path.splitext(fn)
    ext = ext.lower()
    if ext == ".bz2":
        import bz2
        return bz2.open(fn, *args, **kwargs)
    elif ext == ".gz":
        import gzip
        return gzip.open(fn, *args, **kwargs)

    return open(fn, *args, **kwargs)