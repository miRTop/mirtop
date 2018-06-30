from mirtop.gff.body import read_gff_line
import mirtop.libs.logger as mylog
logger = mylog.getLogger(__name__)

def _check_header(header):
    """ Check that header has the minimum mandatory fields
    """
    # Check mandatory fields present:
    num_samples = 0
    matching = []
    mandatory_fields = ["source-ontology","COLDATA","VERSION"]
    for field in mandatory_fields:
            match = ([s for s in header if field in s])
            if len(match) != 0:
                matching.append(match)

    all_present = len(matching)==len(mandatory_fields)
    if all_present:
        for line in header:
            if line.startswith("## COLDATA"):
                num_samples = len(line.strip().split(": ")[1].strip().split(","))


    return [all_present,num_samples]

def _check_line(line,num, num_samples):
    """ Check file for minimum
    """
    fields = read_gff_line(line)
    # Check strand
    if str(line['strand']) not in ["+","-"]:
        print "INCORRECT STRAND"
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
    print("Header obtained")
    all_present, num_samples = _check_header(header)
    if all_present == False:
        raise ValueError("%s doesn't contain all the mandatory fields for the header." % file)
    print("HEADER CHECKED")
    # Check lines

    with open(file) as ch:
        for num, line in enumerate(ch,1):
            if line.startswith("##"):
                next
            else:
                _check_line(line,num,num_samples)

def check_multiple(args):
    for file in args.files:
        _check_file(file)
