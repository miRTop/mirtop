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
    return _guess_database_file(args.gtf)


def _guess_database_file(gff):
    database = None
    with open(gff) as in_handle:
        for line in in_handle:
            if not line.startswith("#"):
                break
            if line.find("miRBase") > -1:
                database = line[line.find("miRBase"):].strip().replace(" ", "")
            elif line.find("MirGeneDB") > -1:
                database = line[line.find("MirGeneDB"):].strip().replace(" ", "")
            elif line.find("microRNAs") > -1:
                database = line.strip().split()[1]
    if not database:
        logger.error("Database not found in --mirna %s. "
                     "Use --database argument to add a custom source." % gff)
        raise ValueError("Database not found in %s header" % gff)
    return database


def get_primary_transcript(database):
    """
    Get the ID to identify the primary transcript in the
      GTF file with the miRNA and precursor coordinates
      to be able to parse BAM files with genomic
      coordinates.
    """
    if database.find("miRBase") > -1:
        return "miRNA_primary_transcript"
    else:
        raise ValueError("Only miRBase is supported for this action.")


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
    if _guess_database_file(gtf).find("miRBase") > -1:
        mapped = read_gtf_to_precursor_mirbase(gtf, format="genomic")
    elif _guess_database_file(gtf).find("MirGeneDB") > -1:
        mapped = read_gtf_to_precursor_mirgenedb(gtf, format="genomic")
    else:
        logger.info("Database different than miRBase or MirGeneDB")
        logger.info("If you get an error when loading,")
        logger.info("report it to https://github.com/miRTop/mirtop/issues")
        mapped = read_gtf_to_precursor_mirbase(gtf)
    return mapped


def read_gtf_chr2mirna(gtf):  # read from read_gtf_to_precursor
    """
    Load GTF file with precursor positions on genome.

    Args:
        *gtf(str)*: file name with GFF miRNA genomic positions and
            header lines.

    Returns:
        *db_mir(dict)*: dictionary with keys being chr and values
            mirna and genomic positions.
    """
    if not gtf:
        return gtf
    if _guess_database_file(gtf).find("miRBase") > -1:
        mapped = read_gtf_to_precursor_mirbase(gtf, format="chrom")
    elif _guess_database_file(gtf).find("MirGeneDB") > -1:
        mapped = read_gtf_to_precursor_mirgenedb(gtf, format="chrom")
    else:
        logger.info("Database different than miRBase or MirGeneDB")
        logger.info("If you get an error when loading,")
        logger.info("report it to https://github.com/miRTop/mirtop/issues")
        mapped = read_gtf_to_precursor_mirbase(gtf)
    return mapped


def read_gtf_chr2mirna2(gtf):  # to remove
    mir2hairpin = read_gtf_to_precursor(gtf)
    if not gtf:
        return gtf
    db = defaultdict(list)
    db_mir = defaultdict(list)
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
            id_dict[idname[0]] = name[0] # MIMA to sps-Y-X-5/3p
            if cols[2] == "miRNA_primary_transcript":
                db[idname[0]] = [chrom, int(start), int(end), strand]
            if cols[2] == "miRNA":
                parent = [n.split("=")[1] for n in cols[-1].split(";")
                          if n.startswith("Derives_from")]
                parent_name = id_dict[parent[0]]
                db_mir[chrom].append([name[0], start, end,
                                      strand,
                                      parent_name,
                                      mir2hairpin[parent_name][name[0]][0]])
                logger.debug("MAP:: mirna:%s" % name[0])
                logger.debug("MAP:: precursor:%s" % parent_name)
                logger.debug("MAP:: precursor pos %s" % mir2hairpin[parent_name][name[0]])
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
    if _guess_database_file(gtf).find("miRBase") > -1:
        mapped = read_gtf_to_precursor_mirbase(gtf)
    elif _guess_database_file(gtf).find("MirGeneDB") > -1:
        mapped = read_gtf_to_precursor_mirgenedb(gtf)
    else:
        logger.info("Database different than miRBase or MirGeneDB")
        logger.info("If you get an error when loading,")
        logger.info("report it to https://github.com/miRTop/mirtop/issues")
        mapped = read_gtf_to_precursor_mirbase(gtf)
    return mapped


def _parse_db_mir_chrom(db_mir, db, id_dict):
    mir2hairpin = _parse_db_mir(db_mir, db, id_dict)
    map_dict = defaultdict(list)
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
        map_dict[db_mir[mir][0]].append([mir[1], db_mir[mir][1], db_mir[mir][2],
                                         db_mir[mir][3],
                                         id_dict[parent],
                                         mir2hairpin[id_dict[parent]][mir[1]][0]])
    return map_dict


def _parse_db_mir(db_mir, db, id_dict):
    map_dict = defaultdict(dict)
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
        logger.debug("MAP:: final:%s %s %s" % (mir[1], start, end))
        map_dict[id_dict[parent]][mir[1]] = [start, end, db_mir[mir][3]]
    return map_dict


def _parse_db_mir_genomic(db_mir, db, id_dict):
    map_dict = defaultdict(dict)
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
        map_dict[mir[1]].update({id_dict[parent]: db_mir[mir]})
    return map_dict


def read_gtf_to_precursor_mirgenedb(gtf, format="precursor"):
    """
    Load GTF file with precursor positions on genome
    Return dict with key being precursor name and
    value a dict of mature miRNA with relative position
    to precursor. For MirGeneDB and similar GFF3 files.

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
    with open(gtf) as in_handle:
        for line in in_handle:
            if line.startswith("#"):
                continue
            cols = line.strip().split("\t")
            idname = [n.split("=")[1] for n in cols[-1].split(";")
                      if n.startswith("ID")][0]
            name = idname
            chrom, start, end, strand = cols[0], cols[3], cols[4], cols[6]
            id_dict[idname] = name
            if cols[2] == "pre_miRNA":
                db[name] = [chrom, int(start), int(end), strand]
            if cols[2] == "miRNA":
                idname_mi = [n.split("=")[1] for n in cols[-1].split(";")
                             if n.startswith("ID")][0]
                parent = "%s_pre" % idname_mi.split("_")[0]
                db_mir[(parent, name)] = [chrom,
                                          int(start), int(end),
                                          strand, parent]
                logger.debug("MAP:: mirna:%s" % name)
                logger.debug("MAP:: pos %s" % db_mir[(parent, name)])
    if format == "genomic":
        return _parse_db_mir_genomic(db_mir, db, id_dict)
    elif format == "chrom":
        return _parse_db_mir_chrom(db_mir, db, id_dict)
    return _parse_db_mir(db_mir, db, id_dict)


def read_gtf_to_precursor_mirbase(gtf, format="precursor"):
    """
    Load GTF file with precursor positions on genome
    Return dict with key being precursor name and
    value a dict of mature miRNA with relative position
    to precursor. For miRBase and similar GFF3 files.

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
    if format == "genomic":
        return _parse_db_mir_genomic(db_mir, db, id_dict)
    elif format == "chrom":
        return _parse_db_mir_chrom(db_mir, db, id_dict)
    return _parse_db_mir(db_mir, db, id_dict)


def liftover_genomic_precursor(read, genome, hairpin, expected=None):
    # example 1 LESS (+) strand to miRNAX (+)1010 (->10)
    # genome miRNAX [start, end, strand] (+)1009
    # hairpin miRNAX [start, end] 9 (1 MORE = 8, 1 LESS = 10)
    # hairpin_start - (genome_start - read_start)
    # 9 - (1009 - 1010) = 9 - (-1) = 10
    if read['strand'] != genome['strand']:
        # Add warning
        return None
    logger.debug("MAPPER::liftover::genomic %s" % genome)
    logger.debug("MAPPER::liftover::hairpin %s" % hairpin)
    logger.debug("MAPPER::liftover::read %s" % read)
    if genome['strand'] == "+":
        lifted = hairpin["start"] - (genome["start"] - read["start"])
    # example 1 MORE: 1008
    # 9 - (1009 - 1008) = 9 - 1 = 8
    # example 1 MORE (-) strand: 1010
    # hairpin_start - (read_start - genome_start)
    # 9 - (1010 - 1009) = 9 - 1 = 8
    # example 1 LESS: 1008
    #  9 - (1008 - 1009) = 9 - (-1) = 10
    elif genome['strand'] == "-":
        lifted = hairpin["start"] - (read["end"] - genome["end"])
    logger.debug("MAPPER::liftover::lifted %s" % lifted)
    if expected and expected != lifted:
        raise ValueError("Bad liftover event.")
    return lifted
