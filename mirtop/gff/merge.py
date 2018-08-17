from collections import defaultdict

from mirtop.gff.body import paste_columns, guess_format
import mirtop.libs.logger as mylog
from mirtop.gff.classgff import feature

logger = mylog.getLogger(__name__)


def merge(dts, samples):
    """
    For dictionary with sample as keys and values as lines
    merge them into one GFF file.

    Args:
        *dts(dict)*: dictionary as >>> {'file': {'mirna': {start: gff_list}}}.
        gff_list has the format as defined in *mirtop.gff.body.read()*.

        *samples(list)*: character list with sample names.

    Returns:
        *merged_lines (nested dicts)*:gff_list has the format as defined in *mirtop.gff.body.read()*.
    """
    logger.debug("MERGE::SAMPLES::given %s" % samples)
    all_data = defaultdict(dict)
    all_lines = defaultdict(list)
    merged_lines = defaultdict(dict)
    for fn in dts:
        for m in dts[fn]:
            for s in dts[fn][m]:
                for hit in dts[fn][m][s]:
                    idu = hit[0]
                    logger.debug("MERGE::SAMPLES::counts %s" % [hit[3], hit[2]])
                    formatted_counts = _format_samples_counts(hit[3], hit[2])
                    logger.debug("MERGE::SAMPLES::fixed %s" % formatted_counts)
                    for sample in formatted_counts:
                        all_data[idu][sample] = formatted_counts[sample] # get the expression of the sample from line
                    all_lines[idu] = hit[4] # get the line
    for idu in all_data:
        line = all_lines[idu]
        expression = _convert_to_string(all_data[idu], samples)
        if _start(line) not in merged_lines[_chrom(line)]:
            merged_lines[_chrom(line)][_start(line)] = []
        merged_lines[_chrom(line)][_start(line)].append([idu, "", "", "",
                                                             _fix(line, expression)])
    return merged_lines


def _format_samples_counts(samples, expression):
    """Return a dictionary of samples counts"""
    if isinstance(samples, list):
        if len(samples) != len(expression):
            raise ValueError("samples %s has different length than expression %s" % (samples,
                                                                                     expression))
    else:
        samples = [samples]
        expression = [expression]
    return dict(zip(samples, expression))


def _fix(line, expression):
    # Need to fix Read attribute since not usefull when multiple sample in a line.
    gff = feature(line)
    attr = gff.attributes
    attr['Expression'] = expression
    return paste_columns(gff, guess_format(line))


def _convert_to_string(d, s):
    v = [str(d[ss]) if ss in d else "0" for ss in s]
    return "%s" % ",".join(v)


def _chrom(string):
    return string.split("\t")[0]


def _start(string):
    return string.split("\t")[3]
