from __future__ import print_function

import datetime
import sys
import os.path as op

import six

from mirtop.mirna.fasta import read_precursor
from mirtop.mirna.mapper import read_gtf_to_precursor, read_gtf_to_mirna
from mirtop.gff.body import read_gff_line
import mirtop.libs.logger as mylog

logger = mylog.getLogger(__name__)


def convert(args):
    for fn in args.files:
        out_file = op.join(args.out, "%s.vcf" % op.splitext(op.basename(fn))[0])
        logger.info("Reading %s" % fn)
        create_vcf(fn, args.hairpin, args.gtf, out_file)
        logger.info("VCF generated %s" % out_file)


def cigar_2_key(cigar, readseq, refseq, pos, var5p, var3p, parent_ini_pos, parent_end_pos, hairpin):
    """
    Args:
        'cigar(str)': CIGAR standard of a compressed alignment representation, this CIGAR omits the '1' integer.
        'readseq(str)': the read sequence
        'refseq(str)': the reference sequence
        'pos(str)': the start current position
        'var5p(int)': extra nucleotides not in the reference miRNA (5p strand)
        'var3p(int)': extra nucleotides not in the reference miRNA (3p strand)
        'parent_ini_pos(int)': the start position of the parent miRNA
        'parent_end_pos(int)': the last position of the parent miRNA
        'hairpin(str)': the string of the hairpin for all the miRNA
    Returns:
        'key_pos(str list)': a list with the positions of the variants.
        'key_var(str list)': a list with the variant keys found.
        'ref(str)': reference base(s).
        'alt(str)': altered base(s).
    """
    key_pos = []
    key_var = []
    ref = []
    alt = []
    n_Mpar = "0"
    n_M = 0
    n_NM = 0
    n_D = 0
    n_I = 0
    # Obtaining the part of the hairpin outside of its parent region (previous)
    if var5p < 0:
        previous_pos = hairpin[parent_ini_pos + var5p:parent_ini_pos]
        refseq = previous_pos + refseq
    else:
        refseq = refseq[var5p:]
    # Obtaining the part of the hairpin outside of its parent region (posterior)
    if var3p >= 0:
        posterior_pos = hairpin[parent_end_pos+1:parent_end_pos+var3p+1]
        refseq = refseq + posterior_pos
    else:
        refseq = refseq[:var3p]
    # Calculating the type of variant and its ref/alt positions.
    logger.debug("VCF::precursor::read %s" % refseq)
    logger.debug("VCF::CIGAR::read %s" % cigar)
    for i in range(len(cigar)):
        if cigar[i].isdigit():
            n_Mpar = n_Mpar + (cigar[i])
        elif cigar[i] in ["A", "T", "C", "G"]:
            pos_rel_ref = n_M + n_D + n_NM + 1
            key_pos.append(pos + pos_rel_ref - 1)
            key_var.append(cigar[i])
            ref.append(refseq[pos_rel_ref - 1])
            alt.append(cigar[i])
            n_NM += 1
        elif cigar[i] == "D":
            if i == 0:
                print("Unexpected error: 'D' in first Cigar pos")
            elif cigar[i-1] == "D":
                ref[-1] = ref[-1] + refseq[n_M + n_NM + n_D]  # Adds new del in the REF column
                key_var[-1] = "D" + str(int(key_var[-1][1:]) + 1)  # Adds '1' to the number of Dels in succession
                n_D += 1
            else:
                pos_rel_ref = n_M + n_NM + n_D
                key_pos.append(pos + pos_rel_ref)
                ref.append(refseq[pos_rel_ref-1:pos_rel_ref+1])
                alt.append(refseq[pos_rel_ref-1])
                key_var.append("D1")
                n_D += 1
        elif cigar[i] == "I":
            if i == 0:
                print("Unexpected error: 'I' in first Cigar pos")
            elif cigar[i-1] == "I":
                alt[-1] = alt[-1] + readseq[n_M + n_NM + n_I]  # Adds the new Insert in the ALT column
                key_var[-1] = key_var[-1] + 'I' + readseq[n_M + n_NM + n_I]  # Adds new Ins in the Key
                n_I += 1
            else:
                pos_rel_read = n_M + n_NM + n_I
                pos_rel_ref = n_M + n_NM + n_D
                key_pos.append(pos + pos_rel_ref)
                alt.append(readseq[pos_rel_read-1:pos_rel_read+1])
                ref.append(refseq[pos_rel_ref-1])
                key_var.append("I" + alt[-1][-1])
                n_I += 1
        elif cigar[i] == "M":
            if (i == 0) | (not cigar[i - 1].isdigit()):
                n_M = n_M + 1
            else:
                n_M = n_M + int(n_Mpar)
                n_Mpar = "0"
        else:
            print("Unexpected error: %s, index: %s" % (cigar[i], i))
    return(key_pos, key_var, ref, alt)


def create_vcf(mirgff3, precursor, gtf, vcffile):
    """
    Args:
        'mirgff3(str)': File with mirGFF3 format that will be converted
        'precursor(str)': Fasta format sequences of all miRNA hairpins
        'gtf(str)': Genome coordinates
        'vcffile': name of the file to be saved
    Returns:
        Nothing is returned, instead, a VCF file is generated
    """
    #Check if the input files exist:
    try:
        gff3_file = open(mirgff3, "r", encoding="utf-8") if six.PY3 else open(mirgff3, "r")
    except IOError:
        print ("Can't read the file", end=mirgff3)
        sys.exit()
    with gff3_file:
        data = gff3_file.read()
        if six.PY2:
            data = data.decode("utf-8-sig").encode("utf-8")

    gff3_data = data.split("\n")
    vcf_file = open(vcffile, "w")

    ver = "v4.3"  # Current VCF version formatting
    vcf_file.write("##fileformat=VCF%s\n" % ver)
    date = datetime.datetime.now().strftime("%Y%m%d")
    vcf_file.write("##fileDate=%s\n" % date)
    source = "\n".join(s for s in gff3_data if "## source-ontology: " in s)[20:]
    line = 0
    sample_names = []
    while gff3_data[line][:2] == "##":
        if gff3_data[line][:19] == "## source-ontology:":
            source = gff3_data[line][20:]
        elif gff3_data[line][:11] == "## COLDATA:":
            sample_names = gff3_data[line][12:].split(",")
        line += 1
    vcf_file.write("##source=%s\n" % source)
    vcf_file.write('##INFO=<ID=NS,Type=Integer,Description="Number of samples"\n')
    vcf_file.write("##FILTER=<ID=REJECT,Description='"'Filter not passed'"'>\n")
    vcf_file.write('##FORMAT=<ID=TRC,Number=1,Type=Integer,Description="Total read count">\n')
    vcf_file.write('##FORMAT=<ID=TSC,Number=1,Type=Integer,Description="Total SNP count">\n')
    vcf_file.write('##FORMAT=<ID=TMC,Number=1,Type=Integer,Description="Total miRNA count">\n')
    vcf_file.write('##FORMAT=<ID=GT,Number=1,Type=Integer,Description="Genotype">\n')
    header = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT"
    # Adds Header
    for s in range(len(sample_names)):
        header = header + "\t" + sample_names[s]
    vcf_file.write(header)

    all_dict = dict()  # initializing an empty dictionary where all info will be added
    key_list = []  # Initializing a list which will contain all the keys of the dictionary
    mirna_dict = dict()  # initializing an empty dictionary where mirna info will be put
    n_SNP = 0
    n_noSNP = 0
    no_var = 0
    hairpins = read_precursor(precursor)
    gff3 = read_gtf_to_precursor(gtf)
    gtf_dic = read_gtf_to_mirna(gtf)
    for line in range(0, len(gff3_data)):
        if not gff3_data[line]:
            continue
        if gff3_data[line][1] == "#":
            continue
        else:   # Parsing the gff3 mirna lecture:
            gff_fields = read_gff_line(gff3_data[line])
            gtf_name = gff_fields['attrb']['Name']
            gtf_parent = gff_fields['attrb']['Parent']
            if gtf_parent not in gff3:
                continue
            if gtf_name not in gff3[gtf_parent]:
                    continue
            parent_ini_pos = gff3[gtf_parent][gtf_name][0]
            parent_end_pos = gff3[gtf_parent][gtf_name][1]
            ref_seq = (hairpins[gtf_parent][parent_ini_pos:parent_end_pos+1])
            vcf_chrom = gtf_dic[gtf_name][gtf_parent][0]
            vcf_pos = int(gff_fields['start']) + int(gtf_dic[gtf_name][gtf_parent][1])
            hairpin = hairpins[gtf_parent]
            variants = gff_fields['attrb']['Variant'].split(",")
            logger.debug("VCF::Variant::%s" % variants)
            #  Obtaining the iso_3p, iso_add3p and iso_5p values:

            var3p = [s for s in variants if 'iso_3p' in s]
            if len(var3p):
                var3p = int(var3p[0][7:])  # Position of iso_3p value
            else:
                var3p = 0

            var_add3p = [s for s in variants if 'iso_add3p' in s]
            if len(var_add3p):
                var_add3p = int(var_add3p[0][10:])  # Position of iso_add3p value
            else:
                var_add3p = 0
            var3p = var3p + var_add3p
            logger.debug("VCF::VAR_3p::%s" % var3p)
            var5p = [s for s in variants if 'iso_5p' in s]
            if len(var5p):
                var5p = int(var5p[0][7:])  # Position of iso_5p value
            else:
                var5p = 0  #
            logger.debug("VCF::VAR_5p::%s" % var5p)
            cigar = gff_fields['attrb']["Cigar"]
            # Obtaining all the variants from the cigar:
            if 1:
                (key_pos, key_var, vcf_ref, vcf_alt) = cigar_2_key(cigar, gff_fields['attrb']['Read'], ref_seq, vcf_pos,
                                                                   var5p, var3p, parent_ini_pos, parent_end_pos, hairpin)

                # Adding the variants to a dictionary and calculating all the fields of a vcf file format:
                if len(key_var) > 0:
                    for s in range(len(key_var)):
                        key_dict = vcf_chrom + '-' + str(key_pos[s]) + '-' + str(key_var[s])
                        raw_counts = gff_fields['attrb']['Expression']
                        raw_counts = [int(i) for i in raw_counts.split(',')]
                        nozero_counts = [int(i > 0) for i in raw_counts]  # counts for every sample if expr != 0.
                        if gtf_name in mirna_dict:  # Adding expression values to same mirnas
                            mirna_dict[gtf_name]['Z'] = [sum(x) for x in zip(mirna_dict[gtf_name]['Z'], raw_counts)]
                        else:
                            mirna_dict[gtf_name] = {}
                            mirna_dict[gtf_name]["Z"] = raw_counts
                        if key_dict in all_dict:
                            if all_dict[key_dict]["Type"] in ["A", "C", "T", "G"]:
                                all_dict[key_dict]['X'] = [sum(x) for x in zip(all_dict[key_dict]['X'], nozero_counts)]
                                all_dict[key_dict]['Y'] = [sum(x) for x in zip(all_dict[key_dict]['Y'], raw_counts)]
                        else:
                            key_list.append(key_dict)
                            all_dict[key_dict] = {}
                            all_dict[key_dict]["Chrom"] = vcf_chrom
                            all_dict[key_dict]["Position"] = key_pos[s]
                            all_dict[key_dict]["mirna"] = gtf_name
                            all_dict[key_dict]["Type"] = key_var[s]
                            if key_var[s][0] in ["A", "C", "T", "G"]:
                                n_SNP += 1
                                all_dict[key_dict]["SNP"] = True
                                all_dict[key_dict]["ID"] = gff_fields['attrb']['Name'] + '-SNP' + str(n_SNP)
                                all_dict[key_dict]['X'] = nozero_counts
                                all_dict[key_dict]['Y'] = raw_counts
                            else:
                                n_noSNP += 1
                                all_dict[key_dict]["SNP"] = False
                                all_dict[key_dict]["ID"] = gff_fields['attrb']['Name'] + '-nonSNP' + str(n_noSNP)
                            all_dict[key_dict]["Ref"] = vcf_ref[s]
                            all_dict[key_dict]["Alt"] = vcf_alt[s]
                            all_dict[key_dict]["Qual"] = "."
                            all_dict[key_dict]["Filter"] = gff_fields['attrb']['Filter']
                            all_dict[key_dict]["Info"] = "NS=" + str(len(sample_names))
            else:
                no_var += 1

    #  Writing the VCF file:
    for s in key_list:
        variant_line = ("\n%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s" %
                       (all_dict[s]["Chrom"], all_dict[s]["Position"], all_dict[s]["ID"],
                        all_dict[s]["Ref"], all_dict[s]["Alt"], all_dict[s]["Qual"],
                        all_dict[s]["Filter"], all_dict[s]["Info"]))
        if all_dict[s]["Type"] in ["A", "T", "C", "G"]:
            format_col = "TRC:TSC:TMC:GT"
            variant_line = variant_line + "\t" + format_col
            samples = ""
            for n in range(len(sample_names)):
                X = all_dict[s]["X"][n]
                Y = all_dict[s]["Y"][n]
                Z = mirna_dict[all_dict[s]["mirna"]]["Z"][n]
                # Calculating the genotype:
                if Y == 0:
                    GT = "0|0"
                elif Z == Y:
                    GT = "1|1"
                else:
                    GT = "1|0"
                samples = samples + "\t" + str(X) + ":" + str(Y) + ":" + str(Z) + ":" + GT
            variant_line = variant_line + samples
        else:
            format_col = ""
            variant_line = variant_line + format_col
        vcf_file.write(variant_line)
    vcf_file.close()
