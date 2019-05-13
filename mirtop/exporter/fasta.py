""" Read GFF files and output FASTA format"""
from __future__ import print_function

import os
import sys

import mirtop.libs.logger as mylog
from mirtop.gff.classgff import feature
from mirtop.mirna.realign import read_id

logger = mylog.getLogger(__name__)


def convert(args):
    """
    Main function to convert from GFF3 to FASTA format.

    Args:
      *args*: supported options for this sub-command.
        See *mirtop.libs.parse.add_subparser_export()*.
    """
    for fn in args.files:
        logger.info("Reading %s" % fn)
        _process(fn, args.out)


def _process(fn, out_dir):
    if out_dir:
        out_fasta = os.path.join(out_dir, "%s.fasta" %
                                 os.path.splitext(os.path.basename(fn))[0])
    outh = sys.stdout if not out_dir else open(out_fasta, 'w')
    with open(fn) as inh:
        for line in inh:
            if line.startswith("#"):
                continue
            gff = feature(line)
            attr = gff.attributes
            read = read_id(attr["UID"])
            print((">{0}\n{1}").format(attr["UID"], read), file=outh)
