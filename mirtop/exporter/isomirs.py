""" Read GFF files and output isomiRs compatible format"""
from __future__ import print_function

import os

import mirtop.libs.logger as mylog
from mirtop.mirna import fasta, mapper
from mirtop.gff.classgff import feature
from mirtop.gff.header import read_samples
from mirtop.mirna.realign import get_mature_sequence, align_from_variants
from mirtop.mirna.realign import read_id, variant_to_5p, \
                                 variant_to_3p, variant_to_add

logger = mylog.getLogger(__name__)


def convert(args):
    """
    Main function to convert from GFF3 to isomiRs Bioc Package.

    Args:
      *args*: supported options for this sub-command.
        See *mirtop.libs.parse.add_subparser_export()*.
    """
    precursors = fasta.read_precursor(args.hairpin, args.sps)
    matures = mapper.read_gtf_to_precursor(args.gtf)
    for fn in args.files:
        logger.info("Reading %s" % fn)
        _read_file(fn, precursors, matures, args.out)


def _read_file(fn, precursors, matures, out_dir):
    samples = read_samples(fn)
    for sample in samples:
        with open(os.path.join(out_dir, "%s.mirna" % sample), 'w') as outh:
            print("\t".join(
                ["seq", "name", "freq", "mir", "start", "end",
                 "mism", "add", "t5", "t3", "s5", "s3", "DB",
                 "precursor", "ambiguity"]), file=outh)
    with open(fn) as inh:
        for line in inh:
            if line.startswith("#"):
                continue
            gff = feature(line)
            cols = gff.columns
            attr = gff.attributes
            read = read_id(attr["UID"])
            t5 = variant_to_5p(precursors[attr["Parent"]],
                               matures[attr["Parent"]][attr["Name"]],
                               attr["Variant"])
            t3 = variant_to_3p(precursors[attr["Parent"]],
                               matures[attr["Parent"]][attr["Name"]],
                               attr["Variant"])
            add = variant_to_add(read,
                                 attr["Variant"])
            mature_sequence = get_mature_sequence(
                precursors[attr["Parent"]],
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
            logger.debug("exporter::isomir::decode %s" % [attr["Variant"],
                                                          t5, t3, add, mm])
            # Error if attr["Read"] doesn't exist
            print(cols)
            line = [read, attr["Read"], "0", attr["Name"],
                    cols['source'], cols['type'],
                    mm, add, t5, t3, "NA", "NA", "miRNA", attr["Parent"], hit]
            for sample, counts in zip(samples, attr["Expression"].split(",")):
                with open(os.path.join(out_dir, "%s.mirna" % sample),
                          'a') as outh:
                    line[2] = counts
                    print("\t".join(line), file=outh)
