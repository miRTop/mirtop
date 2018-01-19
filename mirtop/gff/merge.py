from collections import defaultdict

def merge(dts):
    """
    For dict with samples lines merge into one file
    """
    all_data = defaultdict(dict)
    all_lines = defaultdict()
    merged_lines = defaultdict(dict)
    samples = set()
    for fn in dts:
        for m in dts[fn]:
            for s in dts[fn][m]:
                for hit in dts[fn][m][s]:
                    idu = hit[0]
                    samples.add(hit[3])
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
    cols = line.split("\t")
    attr = cols[8].split(";")
    cols[8] = ";".join([a if a.find("Expression") < 0 else expression for a in attr])
    return "\t".join(cols)

def _convert_to_string(d, s):
    v = [d[ss] if ss in d else 0 for ss in s]
    return " Expression %s" % ",".join(map(str, v))

def _chrom(string):
    return string.split("\t")[0]

def _start(string):
    return string.split("\t")[3]


