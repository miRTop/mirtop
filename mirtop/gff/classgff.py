from mirtop.gff import gff_versions
from collections import OrderedDict


class feature:
    """"Class with alignment information."""

    def __init__(self, line, sep="="):
        # if isinstance(line, basestring) # str in python 3
        if isinstance(line, dict):
            line = self.create_line(line, sep)
        self.line = line
        self.attributes = {}
        self.columns = {}
        self.read_gff_line()
        self.reserved_attributes = gff_versions.ATTRv['1.1']

    def guess_format(self):
        return "=" if self.line.find("Name=") > -1 else " "

    def create_line(fields, sep):
        """
        Create GFF/GTF line from dictionary with fields
        """
        # TODO: need unit test. Check all exists.
        fields["attrb"] = ("Read {seq_name}; UID {idseq}; Name {name};"
                           "Parent {parent}; "
                           "Variant {variant}; Cigar {cigar};"
                           "Expression {counts}; "
                           "Filter {filter}; Hits {hits};").format(**fields)

        line = ("{chrom}\t{database}\t{source}\t{start}\t{end}"
                "\t{score}\t{strand}\t.\t{attrb}").format(**fields)
        return line

    def paste_columns(self, sep=None):
        """
        Create GFF/GTF line from read_gff_line
        """
        sep = self.guess_format()
        attributes = "%s;" % "; ".join(
            "%s%s%s" % (a, sep, self.attributes[a]) for a in self.attributes)
        return "\t".join([self.columns['chrom'], self.columns['source'],
                          self.columns['type'],
                          self.columns['start'], self.columns['end'],
                          self.columns['score'],
                          self.columns['strand'], self.columns['ext'],
                          attributes])

    def read_attributes(self, gff_attrb, sep=" "):
        """
        Read GFF/GTF attributes and return dictionary with fields
        """
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
