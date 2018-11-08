"""Update gff3 files to newest version"""
from __future__ import print_function
import sys
import os

from mirtop.gff.classgff import feature
from mirtop.gff.header import read_version, get_gff_version
from mirtop.mirna.keys import *
from mirtop.mirna.realign import make_id

import mirtop.libs.logger as mylog

logger = mylog.getLogger(__name__)


def convert(args):
    for fn in args.files:
        out_fn = os.path.join(args.out, os.path.basename(fn))
        update_file(fn, out_fn)


def read_uid_10(idu):
    seq = ""
    for i in idu:
        if i == "1" or i == "2":
            return seq[:-int(i)]
        else:
            if i not in CODE2NT:
                logger.error("UID is not valid (%s)" % idu)
                return False
            seq += CODE2NT[i]
    return seq


def to10to11(gff_line):
    gff_line = gff_line.replace("_snp", "_snv")
    gff_line = gff_line.replace("_add", "_add3p")
    features = feature(gff_line)
    if "iso_5p" in features.attributes["Variant"]:
        variants = features.attributes["Variant"].split(",")
        iso_5p = [v.split(":") for v in variants if v.startswith("iso_5p")]
        iso_5p = -1 * int(iso_5p[0][1])
        if iso_5p > 0:
            iso_5p = "+%s" % iso_5p
        variants = ["iso_5p:%s" % iso_5p if v.startswith("iso_5p") else v for v in variants]
        features.attributes["Variant"] = ",".join(variants)
    features.attributes["UID"] = make_id(read_uid_10(features.attributes["UID"]))
    return features.paste_columns()


def update_file(gff_file, new_gff_file):
    """Update file from file version to current version"""
    versions = ["1.0"]
    functions = [to10to11]

    version = read_version(gff_file)
    init = versions.index(version)

    outh = sys.stdout if not new_gff_file else open(new_gff_file, 'w')
    with open(gff_file) as inh:
        for line in inh:
            if line.find("VERSION") > -1:
                print(get_gff_version(), file=outh)
                continue
            if line.startswith("##"):
                print(line.strip(), file=outh)
                continue
            for updates in range(init, len(functions)):
                print(functions[updates](line), file=outh)
