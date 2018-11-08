"""
Compare multiple GFF files to a reference
"""

from __future__ import print_function

import os

from mirtop.gff.classgff import feature
from mirtop.mirna.realign import read_id
import mirtop.libs.logger as mylog

logger = mylog.getLogger(__name__)


def compare(args):
    """
    From a list of GFF files produce comparison with a reference set.

    Args:
        *args(namedtuple)*: arguments parsed from command line with
            *mirtop.libs.parse.add_subparser_compare()*.
            First file will be considered the reference set.

    Returns:
        *(out_file)*: comparison of the GFF files with the reference.
    """
    result = dict()
    reference = read_reference(args.files[0])
    for fn in args.files[1:]:
        if not os.path.exists(fn):
            raise IOError("%s doesn't exist" % fn)
        result[os.path.basename(fn)] = _compare_to_reference(fn, reference)
    if args.out != "tmp_mirtop":
        fn_out = os.path.join(args.out, "summary.txt")
        with open(fn_out, 'w') as outh:
            for fn in result:
                print("sample\tidu\tseq\ttag\tsame_mirna\t%s" % "\t".join(result[fn][0][3].keys()), file=outh)
                for line in result[fn]:
                    read = read_id(line[0])
                    acc = "\t".join([line[3][v] for v in line[3]])
                    print("%s\t%s\t%s\t%s\t%s\t%s" % (fn, line[0], read, line[1], line[2], acc), file=outh)


def read_reference(fn):
    """Read GFF into UID:Variant

    Args:
        *fn (str)*: GFF file.

    Returns:
        *srna (dict)*: dict with >>> {'UID': 'iso_snp:-2,...'}
    """
    srna = dict()
    with open(fn) as inh:
        for line in inh:
            if line.startswith("#"):
                continue
            gff = feature(line)
            attr = gff.attributes
            srna[attr['UID']] = [_simplify(attr['Variant']), attr]
    return srna


def _compare_to_reference(fn, reference):
    same = 0
    diff = list()
    extra = list()
    miss = list()
    results = list()
    seen = 0
    seen_reference = set()
    with open(fn) as inh:
        for line in inh:
            if line.startswith("#"):
                continue
            gff = feature(line)
            attr = gff.attributes
            if attr['UID'] in reference:
                mirna = "Y" if attr['Name'] == reference[attr['UID']][1]['Name'] else attr['Name']
                accuracy = _accuracy(_simplify(attr['Variant']), reference[attr['UID']][0])
                results.append([attr['UID'], "D", mirna, accuracy])
                if _simplify(attr['Variant']) == reference[attr['UID']][0]:
                    same += 1
                else:
                    diff.append("%s | reference: %s" % (line.strip(), reference[attr['UID']][1]))
                seen += 1
                seen_reference.add(attr['UID'])
            else:
                extra.append("%s | extra" % line.strip())
                results.append([attr['UID'], "E", attr['Name'], _accuracy(_simplify(attr['Variant']), "")])
    for uid in reference:
        if uid not in seen_reference:
            results.append([uid, "M", "N", _accuracy("", reference[uid][0])])
            miss.append("| miss %s" % reference[uid][1])
    logger.info("Number of sequences found in reference: %s" % seen)
    logger.info("Number of sequences matches reference: %s" % same)
    logger.info("Number of sequences different than reference: %s" % len(diff))
    logger.info("Number of sequences extra sequences: %s" % len(extra))
    logger.info("Number of sequences missed sequences: %s" % len(miss))
    return results


def _simplify(variant):
    simple = [v.split(":")[0] for v in variant.split(",")]
    return ",".join(simple)


def _get_samples(fn):
    with open(fn) as inh:
        for line in inh:
            if line.startswith("## COLDATA"):
                return line.strip().split(": ")[1].strip().split(",")
    raise ValueError("%s doesn't contain COLDATA header." % fn)


def _accuracy(target, reference):
    """Compare each isomir label in Variant field
       and return a list with values whether:
           FP: no in reference
           FN: no in target
           TP: same values
    """
    logger.debug("COMPARE::ACCURACY::values %s vs %s" % (target, reference))
    accuracy = dict()
    types = ["iso_5p", "iso_3p", "iso_add3p", "iso_snv",
             "iso_snv_seed", "iso_snv_central",
             "iso_snv_central_supp", "iso_snv_central_offset"]
    for t in types:
        if t in reference:
            accuracy[t] = "TP" if t in target else "FN"
        else:
            accuracy[t] = "TN" if t not in target else "FP"
    logger.debug("COMPARE::ACCURACY::%s" % accuracy.keys())
    logger.debug("COMPARE::ACCURACY::%s" % accuracy.values())
    return accuracy
