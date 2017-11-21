from collections import defaultdict

def merge(dts):
    """
    For dict with samples lines merge into one file
    """
    all_data = defaultdict(dict)
    all_lines = defaultdict()
    for fn in dts:
        for m in dts[fn]:
            for s in dts[fn][m]:
                idu = dts[fn][m][0]
                all_data[idu][dts[fn][m][3]] = dts[fn][m][1]
                all_lines[idu] = dts[fn][m][4]
    for idu in all_data:
        expression = _convert_to_string(all_data[idu])
        if _start(all_lines[idu]) not in merged_lines[_chrom(all_lines[idu])]:
            merged_lines[_chrom(all_lines[idu])][_start(all_lines[idu])] = []
        merged_lines[_chrom(all_lines[idu])][_start(all_lines[idu])] = _fix(all_lines[idu], expression)

def _fix(line, expression):
    cols = line.split("\t")
    attr = cols[8].split(";")
    cols[8] = ";".join([a if a.find("Expression") < 0 else expression for a in attr])
    return "\t".join(cols)

def _convert_to_string(d, s):
    v = [d[ss] for ss in s]
    return " Expression %s" % ",".join(map(str, v))

def _chrom(string):
    return string.split("\t")[0]

def _start(string):
    return string.split("\t")[3]


