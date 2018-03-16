""" Read GFF files and output isomiRs compatible format"""

import traceback
import os.path as op
import os
import re
import shutil
from collections import defaultdict

from mirtop.libs import do
from mirtop.libs.utils import file_exists
import mirtop.libs.logger as mylog
from mirtop.mirna import fasta, mapper
from mirtop.mirna.realign import isomir, hits
from mirtop.gff.body import read_attributes
from mirtop.gff.header import read_samples
from mirtop.mirna.realign import get_mature_sequence, align_from_variants, read_id

logger = mylog.getLogger(__name__)

def convert(args):
    samples = []
    database = mapper.guess_database(args.gtf)
    precursors = fasta.read_precursor(args.hairpin, args.sps)
    matures = mapper.read_gtf_to_precursor(args.gtf)
    for fn in args.files:
        logger.info("Reading %s" % fn)
        read_file(fn, precursors, matures, args.out)

def read_file(fn, precursors, matures, out_dir):
    samples = read_samples(fn)
    for sample in samples:
        with open(os.path.join(out_dir, "%s.mirna" % sample), 'w') as outh:
            print >>outh, "\t".join(["seq", "name", "freq", "mir", "start", "end",
                                     "mism", "add", "t5", "t3", "s5", "s3", "DB",
                                     "precursor", "ambiguity"])
    with open(fn) as inh:
        for line in inh:
            if line.startswith("#"):
                continue
            cols = line.strip().split("\t")
            attr = read_attributes(line)
            read = read_id(attr["UID"])
            t5 = _get_5p(precursors[attr["Parent"]],
                         matures[attr["Parent"]][attr["Name"]],
                         attr["Variant"])
            t3 = _get_3p(precursors[attr["Parent"]],
                         matures[attr["Parent"]][attr["Name"]],
                         attr["Variant"])
            add = _get_add(read,
                           attr["Variant"])
            mature_sequence = get_mature_sequence(precursors[attr["Parent"]],
                                                  matures[attr["Parent"]][attr["Name"]])
            mm = align_from_variants(read,
                                     mature_sequence,
                                     attr["Variant"])
            if len(mm) > 1:
                continue
            elif len(mm) == 1:
                mm = "".join(map(str, mm[0]))
            else:
                mm = "0"
            hit = attr["Hits"] if "Hits" in attr else "1"
            logger.debug("exporter::isomir::decode %s" % [attr["Variant"], t5, t3, add, mm])
# Error if attr["Read"] doesn't exist
            line = [read, attr["Read"], "0", attr["Name"], cols[1], cols[2], mm, add, t5, t3, "NA", "NA", "miRNA",  attr["Parent"], hit]
            for sample, counts in zip(samples, attr["Expression"].split(",")):
                with open(os.path.join(out_dir, "%s.mirna" % sample), 'a') as outh:
                    line[2] = counts
                    print >>outh, "\t".join(line)

def _get_5p(hairpin, pos, variant):
    """from a sequence and a start position get the nts
       +/- indicated by variants:
           AAATTTT, 3, -1
           return: T
       pos option is 0-base-index
    """
    pos = pos[0]
    t5 = [v for v in variant.split(",") if v.startswith("iso_5p")]
    if t5:
        t5 = int(t5[0].split(":")[-1][-1])
        if t5 > 0:
            return hairpin[pos - t5:pos]
        elif t5 < 0:
            return hairpin[pos:pos + t5].lower()
    return "0"

def _get_3p(hairpin, pos, variant):
    """from a sequence and a end position get the nts
       +/- indicated by variants:
           AAATTTT, 3, -1
           return: A
       pos option is 0-base-index
    """
    pos = pos[1]
    t3 = [v for v in variant.split(",") if v.startswith("iso_3p")]
    if t3:
        t3 = int(t3[0].split(":")[-1][-1])
        if t3 > 0:
            return hairpin[pos:pos + t3]
        elif t3 < 0:
            return hairpin[pos - t3:pos].lower()
    return "0"

def _get_add(read, variant):
    """from a sequence and a end position get the nts
       +/- indicated by variants:
           AAATTTT, 3, 2
           return: TT
       pos option is 0-base-index
       variant is always positive
    """
    add = [v for v in variant.split(",") if v.startswith("iso_add")]
    if add:
        add = int(add[0].split(":")[-1][-1]) * -1
        return read[add:]
    return "0"

def _get_change(alignment):
    """
    from a 2 string alignment, return the positions
    at which the nucleotides are different in a list:
        AAAAGTTT, AAAACTTT
        return 4CG
    """
    return None
