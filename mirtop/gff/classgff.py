from mirtop.gff import gff_versions
from collections import OrderedDict


class feature:
    """"Class with alignment information."""

    def __init__(self, line):
        self.line = line
        self.attributes = {}
        self.columns = {}
        self.read_gff_line()
        self.reserved_attributes = gff_versions.ATTRv['1.0']

    def guess_format(self):
        return "=" if self.line.find("Name=") > -1 else " "

    def paste_columns(cols, sep=" "):
        """
        Create GFF/GTF line from read_gff_line
        """
        cols['attrb'] = "; ".join(
            "%s%s%s" % (a, sep, cols['attrb'][a]) for a in cols['attrb'])
        return "\t".join([cols['chrom'], cols['source'], cols['type'],
                          cols['start'], cols['end'], cols['score'],
                          cols['strand'], cols['ext'], cols['attrb']])

    def read_attributes(self, gff_attrb, sep=" "):
        sep = self.guess_format()
        gff_dict = OrderedDict()
        for gff_item in gff_attrb.split(";"):
            item_pair = gff_item.strip().split(sep)
            if len(item_pair) > 1:
                gff_dict[item_pair[0].strip()] = item_pair[1].strip()
        self.attributes = gff_dict

    def read_gff_line(self):
        """
        Read GFF/GTF line and return dictionary with fields
        """
        line = self.line
        if line.startswith("#"):
            return line
        cols = line.strip().split("\t")
        if len(cols) != 9:
            raise ValueError("Line has less than 9 elements: %s" % line)
        self.read_attributes(cols[8])
        fields = {'chrom': cols[0],
                  'source': cols[1],
                  'type': cols[2],
                  'start': cols[3],
                  'end': cols[4],
                  'score': cols[5],
                  'strand': cols[6],
                  'ext': cols[7]}
        self.columns = fields
