import traceback
import os.path as op
import os
import re
import shutil
import pandas as pd
import pysam

from mirtop.libs import do
from mirtop.libs.utils import file_exists
import mirtop.libs.logger as mylog

logger = mylog.getLogger(__name__)


def read_gtf(gtf):
    """
    Load GTF file with precursor positions on genome
    """
    if not gtf:
        return gtf
    db = defaultdict(list)
    with open(gtf) as in_handle:
        for line in in_handle:
            if line.startswith("#"):
                continue
            cols = line.strip().split("\t")
            name = [n.split("=")[1] for n in cols[-1].split(";") if n.startswith("Name")]
            chrom, start, end, strand = cols[0], cols[3], cols[4], cols[6]
            if cols[2] == "miRNA_primary_transcript":
                db[name[0]].append([chrom, int(start), int(end), strand])
    return db
