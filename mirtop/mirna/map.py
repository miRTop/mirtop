import traceback
import os.path as op
import os
import re
import shutil
import pandas as pd
import pysam
from collections import defaultdict

from mirtop.libs import do
from mirtop.libs.utils import file_exists
import mirtop.libs.logger as mylog

logger = mylog.getLogger(__name__)


def read_gtf_to_precursor(gtf):
    """
    Load GTF file with precursor positions on genome
    Return dict with key being precursor name and
    value a dict of mature miRNA with relative position
    to precursor.
    """
    if not gtf:
        return gtf
    db = defaultdict(list)
    db_mir = defaultdict(list)
    id_dict = dict()
    map_dict = defaultdict(dict)
    with open(gtf) as in_handle:
        for line in in_handle:
            if line.startswith("#"):
                continue
            cols = line.strip().split("\t")
            name = [n.split("=")[1] for n in cols[-1].split(";") if n.startswith("Name")]
            idname = [n.split("=")[1] for n in cols[-1].split(";") if n.startswith("ID")]
            chrom, start, end, strand = cols[0], cols[3], cols[4], cols[6]
            id_dict[idname[0]] = name[0]
            if cols[2] == "miRNA_primary_transcript":
                db[name[0]] = [chrom, int(start), int(end), strand]
            if cols[2] == "miRNA":
                parent = [n.split("=")[1] for n in cols[-1].split(";") if n.startswith("Derives_from")]
                db_mir[(parent[0], name[0])] = [chrom, int(start), int(end), strand, parent[0]]
                logger.debug(name[0])
                logger.debug(db_mir[(parent[0], name[0])])
    for mir in db_mir:
        parent = db_mir[mir][4]
        precursor = db[id_dict[parent]]
        logger.debug("%s %s %s" % (id_dict[parent], precursor[1], precursor[2]))
        logger.debug("%s %s %s" % (mir, db_mir[mir][1], db_mir[mir][2]))
        if precursor[3] != db_mir[mir][3]:
            logger.warning("%s -> %s" % (id_dict[parent], mir))
            logger.warning("miRNA strand doesn't match with precursor strand: %s - %s" % (db_mir[mir][3], precursor[3]))
            next
        if precursor[0] != db_mir[mir][0]:
            logger.warning("%s -> %s" % (id_dict[parent], mir))
            logger.warning("miRNA chr doesn't match with precursor chr: %s - %s" % (db_mir[mir][0], precursor[0]))
            next
        if precursor[3] == "+":
            start = db_mir[mir][1] - precursor[1]
            end = precursor[2] - db_mir[mir][2]
        if precursor[3] == "-":
            end = db_mir[mir][1] - precursor[1]
            start = precursor[2] - db_mir[mir][2]
        db_mir[mir][1] = start
        db_mir[mir][2] = end
        logger.debug("%s %s %s" % (mir[1], start, end))
        map_dict[id_dict[parent]][mir[1]] = db_mir[mir]
    return map_dict
