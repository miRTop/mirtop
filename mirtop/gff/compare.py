"""
Compare multiple GFF files to a reference
"""

import os
import pandas as pd

from mirtop.gff import header
import mirtop.libs.logger as mylog
logger = mylog.getLogger(__name__)

# Add check first 

def compare(args):
    """
    From a list of files produce stats
    """
    out = list()
    reference = read_reference(args.files[0])
    for fn in args.files[1:]:
        if not os.path.exists(fn):
            raise IOError("%s doesn't exist" %s)
        _compare_to_reference(fn, reference)

def read_reference(fn):
    """Read GFF into UID:Variant key:value dict"""
    srna = dict()
    with open(fn) as inh:
        for line in inh:
            if line.startswith("#"):
                continue
            cols = line.strip().split("\t")
            attr_v = [v.strip().split(" ")[1] for v in cols[8].strip().split(";")[:-1]]
            attr_k = [v.strip().split(" ")[0] for v in cols[8].strip().split(";")[:-1]]
            attr = dict(zip(attr_k, attr_v))
            srna[attr['UID']] = attr['Variant']
    return srna

def _compare_to_reference(fn, reference):
    same = 0
    diff = 0
    new = 0
    with open(fn) as inh:
        for line in inh:
            if line.startswith("#"):
                continue
            cols = line.strip().split("\t")
            cols = line.strip().split("\t")
            attr_v = [v.strip().split(" ")[1] for v in cols[8].strip().split(";")[:-1]]
            attr_k = [v.strip().split(" ")[0] for v in cols[8].strip().split(";")[:-1]]
            attr = dict(zip(attr_k, attr_v))
            if attr['UID'] in reference:
                if attr['Variant'] == reference[attr['UID']]:
                    same += 1
                else:
                    diff +=1
                new += 1
    print "Number of sequences found in reference: %s" % new
    print "Number of sequences matches reference: %s" % same
    print "Number of sequences different than reference: %s" % diff

def _get_samples(fn):
    with open(fn) as inh:
        for line in inh:
            if line.startswith("## COLDATA"):
                return line.strip().split(": ")[1].strip().split(",")
    raise ValueError("%s doesn't contain COLDATA header." % fn)

