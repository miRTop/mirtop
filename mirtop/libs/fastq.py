import os
from collections import Counter
from itertools import product
import gzip


def open_fastq(in_file):
    """ open a fastq file, using gzip if it is gzipped
        (from bcbio package)

    Args:
        *in_file(str)*: file name, including full path.

    Returns:
        *(file)*:File instance to go into for loop.
    """
    _, ext = os.path.splitext(in_file)
    if ext == ".gz":
        return gzip.open(in_file, 'rb')
    if ext in [".fastq", ".fq", ".fasta", ".fa"]:
        return open(in_file, 'r')
    return ValueError("File needs to be fastq|fasta|fq|fa [.gz]")


def is_fastq(in_file):
    """Check whether a file has the expected extension to be a fastq file
        txt, fq, fast + gzip/gz (from bcbio package)
    Args:
        *in_file(str)*: file name, including full path.
    
    Returns:
        *(boolean)*: True or False
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
        (from bcbio package).
    Args:
        *in_file(str)*: file name, including full path.
    
    Returns:
        *base, ext(character list)*: basename of the file and 
            extension of the file.
     
    """
    base, ext = os.path.splitext(fn)
    if ext in [".gz", ".bz2", ".zip"]:
        base, ext2 = os.path.splitext(base)
        ext = ext2 + ext
    return base, ext

