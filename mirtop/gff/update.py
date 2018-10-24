"""Update gff3 files to newest version"""
from __future__ import print_function

from mirtop.gff.classgff import feature
from mirtop.gff.header import read_version, get_gff_version
from mirtop.mirna.keys import *
from mirtop.mirna.realign import make_id

import mirtop.libs.logger as mylog

logger = mylog.getLogger(__name__)


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
    features = feature(gff_line)
    if "iso_5p" in features.attributes:
        features.attributes["iso_5p"] = -1 * int(features.attributes["iso_5p"])
    features.attributes["UID"] = make_id(read_uid_10(features.attributes["UID"]))
    return feature.paste_columns()


versions = ["1.0"]
functions = [to10to11]


def update_file(gff_file, new_gff_file):
    """Update file from file version to current version"""
    version = read_version(gff_file)
    # find position of that version in variable versions
    # apply all the function from that postion and on in functions
    with open(new_gff_file, 'w') as outh:
        with open(gff_file) as inh:
            for line in inh:
                if line.find("VERSION") > -1:
                    print(get_gff_version(), file=outh)
                if line.startswith("##"):
                    print(line.strip(), file=outh)
                # for idx in range(init:len(functions)):
                #   functions[idx](gff_line)
