from collections import defaultdict

from mirtop.gff.body import read_gff_line, paste_columns, guess_format
import mirtop.libs.logger as mylog
logger = mylog.getLogger(__name__)

def merge(dts, samples):
    """
    For dict with samples lines merge into one file
    """
    logger.debug("MERGE::SAMPLES::given %s" % samples)
    all_data = defaultdict(dict)
    all_lines = defaultdict()
    merged_lines = defaultdict(dict)
    for fn in dts:
        for m in dts[fn]:
            for s in dts[fn][m]:
                for hit in dts[fn][m][s]:
                    idu = hit[0]
                    logger.debug("MERGE::SAMPLES::counts %s" % [hit[3], hit[2]])
                    formatted_counts = _format_samples_counts(hit[3], hit[2])
                    logger.debug("MERGE::SAMPLES::fixed %s" % formatted_counts)
                    for s in formatted_counts:
                        all_data[idu][s] = formatted_counts[s] # get the expression of the sample from line
                    all_lines[idu] = hit[4] # get the line
    for idu in all_data:
        expression = _convert_to_string(all_data[idu], samples)
        if _start(all_lines[idu]) not in merged_lines[_chrom(all_lines[idu])]:
            merged_lines[_chrom(all_lines[idu])][_start(all_lines[idu])] = []
        merged_lines[_chrom(all_lines[idu])][_start(all_lines[idu])].append([idu, "", "", "", _fix(all_lines[idu], expression)])
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
    cols = read_gff_line(line)
    cols['attrb']['Expression'] = expression
    return paste_columns(cols, guess_format(line))

def _convert_to_string(d, s):
    v = [str(d[ss]) if ss in d else "0" for ss in s]
    return "%s" % ",".join(v)

def _chrom(string):
    return string.split("\t")[0]

def _start(string):
    return string.split("\t")[3]


