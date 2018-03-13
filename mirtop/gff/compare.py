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
    result = dict()
    reference = read_reference(args.files[0])
    for fn in args.files[1:]:
        if not os.path.exists(fn):
            raise IOError("%s doesn't exist" %s)
        result[os.path.basename(fn)] = _compare_to_reference(fn, reference)
    if args.out != "tmp_mirtop":
        fn_out = os.path.join(args.out, "summary.txt")
        with open(fn_out, 'w') as outh:
            for fn in result:
                for line in result[fn]['different']:
                    print >>outh, "%s\t%s" % (fn, line)
                for line in result[fn]['extra']:
                    print >>outh, "%s\t%s" % (fn, line)
                for line in result[fn]['miss']:
                    print >>outh, "%s\t%s" % (fn, line)

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
            srna[attr['UID']] = [_simplify(attr['Variant']), cols[8]]
    return srna

def _compare_to_reference(fn, reference):
    same = 0
    diff = list()
    extra = list()
    miss = list()
    seen = 0
    seen_reference = set()
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
                if _simplify(attr['Variant']) == reference[attr['UID']][0]:
                    same += 1
                else:
                    diff.append("%s | reference: %s" % (line.strip(), reference[attr['UID']][1]))
                seen += 1
                seen_reference.add(attr['UID'])
            else:
                extra.append("%s | extra" % line.strip())
    for uid in reference:
        if uid not in seen_reference:
            miss.append("| miss %s" %  reference[uid][1])
    logger.info("Number of sequences found in reference: %s" % seen)
    logger.info("Number of sequences matches reference: %s" % same)
    logger.info("Number of sequences different than reference: %s" % len(diff))
    logger.info("Number of sequences extra sequences: %s" % len(extra))
    logger.info("Number of sequences missed sequences: %s" % len(miss))
    return {'different': diff, 'extra': extra, 'miss': miss}

def _simplify(variant):
    simple = [v.split(":")[0] for v in variant.split(",")]
    return ",".join(simple)

def _get_samples(fn):
    with open(fn) as inh:
        for line in inh:
            if line.startswith("## COLDATA"):
                return line.strip().split(": ")[1].strip().split(",")
    raise ValueError("%s doesn't contain COLDATA header." % fn)

