from mirtop.gff.classgff import feature
import mirtop.libs.logger as mylog
from mirtop.gff import gff_versions as version
from mirtop.mirna.realign import read_id

logger = mylog.getLogger(__name__)


def _check_header(header):
    """ Check that header has the minimum mandatory fields
    """
    # Check mandatory fields present:
    num_samples = 0
    matching = []
    mandatory_fields = ["source-ontology", "COLDATA", "VERSION"]
    for field in mandatory_fields:
            match = ([s for s in header if field in s])
            if len(match) != 0:
                matching.append(match)

    all_present = len(matching) == len(mandatory_fields)
    if all_present:
        for line in header:
            if line.startswith("## COLDATA"):
                samples = line.strip().split(": ")[1].strip().split(",")
                num_samples = len(samples)
    return [all_present, num_samples]


def _check_line(line, num, num_samples):
    """ Check file for minimum
    """
    gff = feature(line)
    fields = gff.columns
    attr = gff.attributes

    # Check seqID
    if not fields['chrom']:
        logger.error('MISSING seqID in line %s' % (num))

    # Check source
    source = (fields['source']).lower()
    valid_source = False
    valid_sources = ["mirBase", "mirgeneDB"]
    if (any(s.lower() in source for s in valid_sources)):
        valid_source = True

    if valid_source is False:
        logger.error('INCORRECT SOURCE in line %s' % (num))

    # Check type
    type = fields['type']

    source = (fields['source']).lower()
    valid_type = False
    if type in ["ref_miRNA", "isomiR"]:
        valid_type = True

    if valid_type is False:
        logger.error('INCORRECT TYPE in line %s' % (num))

    # Check start/end
    if not fields['start']:
        logger.error('MISSING START value in line %s' % (num))

    if not fields['end']:
        logger.error('MISSING END value in line %s' % (num))

    # Check strand
    if str(fields['strand']) not in ["+", "-"]:
        logger.error('INCORRECT STRAND in line %s' % (num))

    # Check UID
    if 'UID' not in attr:
        logger.error('UID not found in line %s' % (num))
    else:
        if not read_id(attr['UID']):
            logger.error('UID is not in a correct format in line %s. '
                         'Use mirtop gff to fix this or open an issue.' % num)

    # Check attribute-variant
    variant = (attr['Variant']).lower()
    valid_variant = False
    valid_variants = version.GFFv[version.current]
    if (any(s.lower() in variant for s in valid_variants)):
        valid_variant = True

    if valid_variant is False:
        logger.error('INCORRECT VARIANT type in line %s' % (num))

    # Check attribute-expression

    expression = attr['Expression'].strip().split(",")
    expression = list(filter(None, expression))
    if len(expression) != num_samples:
        logger.error('INCORRECT number of EXPRESSION VALUES \
        in line %s' % (num))


def _check_file(file):
    """
    """
    # Get header to check.
    header = []
    with open(file) as ch:
        for line in ch:
            if line.startswith("##"):
                header.append(line)
            else:
                break
    all_present, num_samples = _check_header(header)
    if all_present is False:
        logger.error("%s doesn't contain all \
        the mandatory fields for the header." % file)
    logger.info("HEADER CHECKED")
    # Check lines

    with open(file) as ch:
        for num, line in enumerate(ch, 1):
            if line.startswith("##"):
                next
            else:
                _check_line(line, num, num_samples)


def check_multiple(args):
    for file in args.files:
        _check_file(file)
        logger.info('%s checked' % file)
