"""Read database information"""

from collections import defaultdict

import mirtop.libs.logger as mylog

logger = mylog.getLogger(__name__)


def guess_database(args):
    """
    Guess database name from GFF file.

    Args:
        *gtf(str)*: file name with GFF miRNA genomic positions and
            header lines.

    Returns:
        *database(str)*: name of the database

    TODO: this needs to be generic to other databases.
    """
    database = None
    with open(args.gtf) as in_handle:
        for line in in_handle:
            if not line.startswith("#"):
                break
            if line.find("miRBase") > -1:
                database = line[line.find("miRBase"):].strip().replace(" ", "")
            elif line.find("microRNAs") > -1:
                database = line.strip().split()[1]
    if not database:
        logger.error("Database not found in --mirna %s. "
                     "Use --database argument to add a custom source." % args.gtf)
        raise ValueError("Database not found in %s header" % args.gtf)
    return database


def read_gtf_to_mirna(gtf):
    """
    Load GTF file with precursor positions on genome.

    Args:
        *gtf(str)*: file name with GFF miRNA genomic positions and
            header lines.

    Returns:
        *db_mir(dict)*: dictionary with keys being mirnas and values
            genomic positions.
    """
    if not gtf:
        return gtf
    db = defaultdict(list)
    db_mir = defaultdict(dict)
    id_dict = dict()
    with open(gtf) as in_handle:
        for line in in_handle:
            if line.startswith("#"):
                continue
            cols = line.strip().split("\t")
            logger.debug("MAP:: line:%s" % cols)
            name = [n.split("=")[1] for n in cols[-1].split(";")
                    if n.startswith("Name")]
            idname = [n.split("=")[1] for n in cols[-1].split(";")
                      if n.startswith("ID")]
            chrom, start, end, strand = cols[0], cols[3], cols[4], cols[6]
            logger.debug("MAP:: idname:%s" % idname)
            logger.debug("MAP:: name:%s" % name)
            id_dict[idname[0]] = name[0]
            if cols[2] == "miRNA_primary_transcript":
                db[idname[0]] = [chrom, int(start), int(end), strand]
            if cols[2] == "miRNA":
                parent = [n.split("=")[1] for n in cols[-1].split(";")
                          if n.startswith("Derives_from")]
                db_mir[name[0]].update({id_dict[parent[0]]: db[parent[0]]})
                logger.debug("MAP:: mirna:%s" % name[0])
                logger.debug("MAP:: precursor:%s" % id_dict[parent[0]])
                logger.debug("MAP:: precursor pos %s" % db[parent[0]])
    return db_mir


def read_gtf_to_precursor(gtf):
    """
    Load GTF file with precursor positions on genome
    Return dict with key being precursor name and
    value a dict of mature miRNA with relative position
    to precursor.

    Args:
        *gtf(str)*: file name with GFF miRNA genomic positions and
            header lines.

    Returns:
        *map_dict(dict)*:

        >>> {'parent': {mirna: [start, end]}}
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
            name = [n.split("=")[1] for n in cols[-1].split(";")
                    if n.startswith("Name")]
            idname = [n.split("=")[1] for n in cols[-1].split(";")
                      if n.startswith("ID")]
            chrom, start, end, strand = cols[0], cols[3], cols[4], cols[6]
            id_dict[idname[0]] = name[0]
            if cols[2] == "miRNA_primary_transcript":
                db[name[0]] = [chrom, int(start), int(end), strand]
            if cols[2] == "miRNA":
                parent = [n.split("=")[1] for n in cols[-1].split(";")
                          if n.startswith("Derives_from")]
                db_mir[(parent[0], name[0])] = [chrom,
                                                int(start), int(end),
                                                strand, parent[0]]
                logger.debug("MAP:: mirna:%s" % name[0])
                logger.debug("MAP:: pos %s" % db_mir[(parent[0], name[0])])
    for mir in db_mir:
        parent = db_mir[mir][4]
        precursor = db[id_dict[parent]]
        logger.debug("MAP::%s %s %s" % (
            id_dict[parent], precursor[1], precursor[2]))
        logger.debug("MAP::%s %s %s" % (mir, db_mir[mir][1], db_mir[mir][2]))
        if precursor[3] != db_mir[mir][3]:
            logger.warning("%s -> %s" % (id_dict[parent], mir))
            logger.warning("miRNA strand doesn't match with precursor strand:"
                           " %s - %s" % (db_mir[mir][3], precursor[3]))
            next
        if precursor[0] != db_mir[mir][0]:
            logger.warning("%s -> %s" % (id_dict[parent], mir))
            logger.warning("miRNA chr doesn't match with precursor chr:"
                           " %s - %s" % (db_mir[mir][0], precursor[0]))
            next
        if precursor[3] == "+":
            start = db_mir[mir][1] - precursor[1]
            end = db_mir[mir][2] - precursor[1]
        if precursor[3] == "-":
            end = precursor[2] - db_mir[mir][1]
            start = precursor[2] - db_mir[mir][2]
        db_mir[mir][1] = start
        db_mir[mir][2] = end
        logger.debug("MAP:: final:%s %s %s" % (mir[1], start, end))
        map_dict[id_dict[parent]][mir[1]] = db_mir[mir][1:3]
    return map_dict
