"""Helpers to work with fastq files"""

import os
from itertools import product
import gzip


def open_fastq(in_file):
    """ open a fastq file, using gzip if it is gzipped
        (from bcbio package)

    Args:
        *in_file(str)*: file name.

    Returns:
        *(File)*: file handler.
    """
    _, ext = os.path.splitext(in_file)
    if ext == ".gz":
        return gzip.open(in_file, 'rb')
    if ext in [".fastq", ".fq", ".fasta", ".fa"]:
        return open(in_file, 'r')
    return ValueError("File needs to be fastq|fasta|fq|fa [.gz]")


def is_fastq(in_file):
    """Check whether file is fastq accepting
        txt, fq and fastq extensions understanding
        compression with gzip: .gzip and .gz
        (copy from bcbio)
    Args:
        *in_file(str)*: file name.

    Returns:
        *(boolean)*: Yes or Not.
    """
    fastq_ends = [".txt", ".fq", ".fastq"]
    zip_ends = [".gzip", ".gz"]
    base, first_ext = os.path.splitext(in_file)
    second_ext = os.path.splitext(base)[1]
    if first_ext in fastq_ends:
        return True
    elif (second_ext, first_ext) in product(fastq_ends, zip_ends):
        return True
    else:
        return False


def splitext_plus(fn):
    """Split on file extensions, allowing for zipped extensions.
        (copy from bcbio)

    Args:
        *fn(str)*: file name.

    Returns:
        *base, ext(str, str)*: basename and extesion.
    """
    base, ext = os.path.splitext(fn)
    if ext in [".gz", ".bz2", ".zip"]:
        base, ext2 = os.path.splitext(base)
        ext = ext2 + ext
    return base, ext
