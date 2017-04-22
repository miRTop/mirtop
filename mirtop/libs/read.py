"""functions for explore tool"""
from __future__ import print_function
from collections import defaultdict

from mirtop.libs.do import find_cmd, run

import pysam
import pybedtools

def get_fasta(bed_file, ref, out_fa):
    """Run bedtools to get fasta from bed file"""
    cmd = "bedtools getfasta -s -fi {ref} -bed {bed_file} -fo {out_fa}"
    run(cmd.format(**locals()))
