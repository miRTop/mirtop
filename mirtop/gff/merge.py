from collections import defaultdict

from mirtop.gff.body import read_gff_line, paste_columns, guess_format

def merge(dts, samples):
    """
    For dict with samples lines merge into one file
    """
    all_data = defaultdict(dict)
    all_lines = defaultdict()
    merged_lines = defaultdict(dict)
    for fn in dts:
        for m in dts[fn]:
            for s in dts[fn][m]:
                for hit in dts[fn][m][s]:
                    idu = hit[0]
                    all_data[idu][hit[3]] = hit[2] # get the expression of the sample from line
                    all_lines[idu] = hit[4] # get the line
    for idu in all_data:
        expression = _convert_to_string(all_data[idu], samples)
        if _start(all_lines[idu]) not in merged_lines[_chrom(all_lines[idu])]:
            merged_lines[_chrom(all_lines[idu])][_start(all_lines[idu])] = []
        merged_lines[_chrom(all_lines[idu])][_start(all_lines[idu])].append([idu, "", "", "", _fix(all_lines[idu], expression)])
    return merged_lines

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


