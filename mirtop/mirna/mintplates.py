#!/usr/bin/env python
"""
MINTplates.py: Terms and conditions located in README_TermsOfUse_MINTplates.txt

This code was created by Venetia Pliatsika, Isidore Rigoutsos, Jeffery Ma, Phillipe Loher
It can be used to create/encode molecular "license-plates" from sequences and to also decode the "license-plates"
back to sequences.  While initially created for tRFs (tRNA fragments), this tool can be used for
any genomic sequences including but not limited to:  tRFs, isomiRs, reference miRNA, etc.
For more information on "license-plates", visit https://cm.jefferson.edu/MINTbase and
refer to publications https://www.ncbi.nlm.nih.gov/pubmed/27153631/ and https://www.ncbi.nlm.nih.gov/pubmed/28220888/.
Contact us at: https://cm.jefferson.edu/contact-us/
"""

import sys
import re
import argparse

encode_hash = {'AAAAA': 'BB', 'AAAAC': 'BD', 'AAAAG': 'B0', 'AAAAT': 'BE', 'AAACA': 'BF', 'AAACC': 'B1',
               'AAACG': 'BH', 'AAACT': 'BI', 'AAAGA': 'B2', 'AAAGC': 'BJ', 'AAAGG': 'BK', 'AAAGT': 'B3',
               'AAATA': 'BL', 'AAATC': 'BM', 'AAATG': 'B4', 'AAATT': 'BN', 'AACAA': 'BO', 'AACAC': 'B5',
               'AACAG': 'BP', 'AACAT': 'BQ', 'AACCA': 'B6', 'AACCC': 'BR', 'AACCG': 'BS', 'AACCT': 'B7',
               'AACGA': 'BU', 'AACGC': 'BV', 'AACGG': 'B8', 'AACGT': 'BW', 'AACTA': 'BX', 'AACTC': 'B9',
               'AACTG': 'BY', 'AACTT': 'BZ', 'AAGAA': 'DB', 'AAGAC': 'DD', 'AAGAG': 'D0', 'AAGAT': 'DE',
               'AAGCA': 'DF', 'AAGCC': 'D1', 'AAGCG': 'DH', 'AAGCT': 'DI', 'AAGGA': 'D2', 'AAGGC': 'DJ',
               'AAGGG': 'DK', 'AAGGT': 'D3', 'AAGTA': 'DL', 'AAGTC': 'DM', 'AAGTG': 'D4', 'AAGTT': 'DN',
               'AATAA': 'DO', 'AATAC': 'D5', 'AATAG': 'DP', 'AATAT': 'DQ', 'AATCA': 'D6', 'AATCC': 'DR',
               'AATCG': 'DS', 'AATCT': 'D7', 'AATGA': 'DU', 'AATGC': 'DV', 'AATGG': 'D8', 'AATGT': 'DW',
               'AATTA': 'DX', 'AATTC': 'D9', 'AATTG': 'DY', 'AATTT': 'DZ', 'ACAAA': '0B', 'ACAAC': '0D',
               'ACAAG': '00', 'ACAAT': '0E', 'ACACA': '0F', 'ACACC': '01', 'ACACG': '0H', 'ACACT': '0I',
               'ACAGA': '02', 'ACAGC': '0J', 'ACAGG': '0K', 'ACAGT': '03', 'ACATA': '0L', 'ACATC': '0M',
               'ACATG': '04', 'ACATT': '0N', 'ACCAA': '0O', 'ACCAC': '05', 'ACCAG': '0P', 'ACCAT': '0Q',
               'ACCCA': '06', 'ACCCC': '0R', 'ACCCG': '0S', 'ACCCT': '07', 'ACCGA': '0U', 'ACCGC': '0V',
               'ACCGG': '08', 'ACCGT': '0W', 'ACCTA': '0X', 'ACCTC': '09', 'ACCTG': '0Y', 'ACCTT': '0Z',
               'ACGAA': 'EB', 'ACGAC': 'ED', 'ACGAG': 'E0', 'ACGAT': 'EE', 'ACGCA': 'EF', 'ACGCC': 'E1',
               'ACGCG': 'EH', 'ACGCT': 'EI', 'ACGGA': 'E2', 'ACGGC': 'EJ', 'ACGGG': 'EK', 'ACGGT': 'E3',
               'ACGTA': 'EL', 'ACGTC': 'EM', 'ACGTG': 'E4', 'ACGTT': 'EN', 'ACTAA': 'EO', 'ACTAC': 'E5',
               'ACTAG': 'EP', 'ACTAT': 'EQ', 'ACTCA': 'E6', 'ACTCC': 'ER', 'ACTCG': 'ES', 'ACTCT': 'E7',
               'ACTGA': 'EU', 'ACTGC': 'EV', 'ACTGG': 'E8', 'ACTGT': 'EW', 'ACTTA': 'EX', 'ACTTC': 'E9',
               'ACTTG': 'EY', 'ACTTT': 'EZ', 'AGAAA': 'FB', 'AGAAC': 'FD', 'AGAAG': 'F0', 'AGAAT': 'FE',
               'AGACA': 'FF', 'AGACC': 'F1', 'AGACG': 'FH', 'AGACT': 'FI', 'AGAGA': 'F2', 'AGAGC': 'FJ',
               'AGAGG': 'FK', 'AGAGT': 'F3', 'AGATA': 'FL', 'AGATC': 'FM', 'AGATG': 'F4', 'AGATT': 'FN',
               'AGCAA': 'FO', 'AGCAC': 'F5', 'AGCAG': 'FP', 'AGCAT': 'FQ', 'AGCCA': 'F6', 'AGCCC': 'FR',
               'AGCCG': 'FS', 'AGCCT': 'F7', 'AGCGA': 'FU', 'AGCGC': 'FV', 'AGCGG': 'F8', 'AGCGT': 'FW',
               'AGCTA': 'FX', 'AGCTC': 'F9', 'AGCTG': 'FY', 'AGCTT': 'FZ', 'AGGAA': '1B', 'AGGAC': '1D',
               'AGGAG': '10', 'AGGAT': '1E', 'AGGCA': '1F', 'AGGCC': '11', 'AGGCG': '1H', 'AGGCT': '1I',
               'AGGGA': '12', 'AGGGC': '1J', 'AGGGG': '1K', 'AGGGT': '13', 'AGGTA': '1L', 'AGGTC': '1M',
               'AGGTG': '14', 'AGGTT': '1N', 'AGTAA': '1O', 'AGTAC': '15', 'AGTAG': '1P', 'AGTAT': '1Q',
               'AGTCA': '16', 'AGTCC': '1R', 'AGTCG': '1S', 'AGTCT': '17', 'AGTGA': '1U', 'AGTGC': '1V',
               'AGTGG': '18', 'AGTGT': '1W', 'AGTTA': '1X', 'AGTTC': '19', 'AGTTG': '1Y', 'AGTTT': '1Z',
               'ATAAA': 'HB', 'ATAAC': 'HD', 'ATAAG': 'H0', 'ATAAT': 'HE', 'ATACA': 'HF', 'ATACC': 'H1',
               'ATACG': 'HH', 'ATACT': 'HI', 'ATAGA': 'H2', 'ATAGC': 'HJ', 'ATAGG': 'HK', 'ATAGT': 'H3',
               'ATATA': 'HL', 'ATATC': 'HM', 'ATATG': 'H4', 'ATATT': 'HN', 'ATCAA': 'HO', 'ATCAC': 'H5',
               'ATCAG': 'HP', 'ATCAT': 'HQ', 'ATCCA': 'H6', 'ATCCC': 'HR', 'ATCCG': 'HS', 'ATCCT': 'H7',
               'ATCGA': 'HU', 'ATCGC': 'HV', 'ATCGG': 'H8', 'ATCGT': 'HW', 'ATCTA': 'HX', 'ATCTC': 'H9',
               'ATCTG': 'HY', 'ATCTT': 'HZ', 'ATGAA': 'IB', 'ATGAC': 'ID', 'ATGAG': 'I0', 'ATGAT': 'IE',
               'ATGCA': 'IF', 'ATGCC': 'I1', 'ATGCG': 'IH', 'ATGCT': 'II', 'ATGGA': 'I2', 'ATGGC': 'IJ',
               'ATGGG': 'IK', 'ATGGT': 'I3', 'ATGTA': 'IL', 'ATGTC': 'IM', 'ATGTG': 'I4', 'ATGTT': 'IN',
               'ATTAA': 'IO', 'ATTAC': 'I5', 'ATTAG': 'IP', 'ATTAT': 'IQ', 'ATTCA': 'I6', 'ATTCC': 'IR',
               'ATTCG': 'IS', 'ATTCT': 'I7', 'ATTGA': 'IU', 'ATTGC': 'IV', 'ATTGG': 'I8', 'ATTGT': 'IW',
               'ATTTA': 'IX', 'ATTTC': 'I9', 'ATTTG': 'IY', 'ATTTT': 'IZ', 'CAAAA': '2B', 'CAAAC': '2D',
               'CAAAG': '20', 'CAAAT': '2E', 'CAACA': '2F', 'CAACC': '21', 'CAACG': '2H', 'CAACT': '2I',
               'CAAGA': '22', 'CAAGC': '2J', 'CAAGG': '2K', 'CAAGT': '23', 'CAATA': '2L', 'CAATC': '2M',
               'CAATG': '24', 'CAATT': '2N', 'CACAA': '2O', 'CACAC': '25', 'CACAG': '2P', 'CACAT': '2Q',
               'CACCA': '26', 'CACCC': '2R', 'CACCG': '2S', 'CACCT': '27', 'CACGA': '2U', 'CACGC': '2V',
               'CACGG': '28', 'CACGT': '2W', 'CACTA': '2X', 'CACTC': '29', 'CACTG': '2Y', 'CACTT': '2Z',
               'CAGAA': 'JB', 'CAGAC': 'JD', 'CAGAG': 'J0', 'CAGAT': 'JE', 'CAGCA': 'JF', 'CAGCC': 'J1',
               'CAGCG': 'JH', 'CAGCT': 'JI', 'CAGGA': 'J2', 'CAGGC': 'JJ', 'CAGGG': 'JK', 'CAGGT': 'J3',
               'CAGTA': 'JL', 'CAGTC': 'JM', 'CAGTG': 'J4', 'CAGTT': 'JN', 'CATAA': 'JO', 'CATAC': 'J5',
               'CATAG': 'JP', 'CATAT': 'JQ', 'CATCA': 'J6', 'CATCC': 'JR', 'CATCG': 'JS', 'CATCT': 'J7',
               'CATGA': 'JU', 'CATGC': 'JV', 'CATGG': 'J8', 'CATGT': 'JW', 'CATTA': 'JX', 'CATTC': 'J9',
               'CATTG': 'JY', 'CATTT': 'JZ', 'CCAAA': 'KB', 'CCAAC': 'KD', 'CCAAG': 'K0', 'CCAAT': 'KE',
               'CCACA': 'KF', 'CCACC': 'K1', 'CCACG': 'KH', 'CCACT': 'KI', 'CCAGA': 'K2', 'CCAGC': 'KJ',
               'CCAGG': 'KK', 'CCAGT': 'K3', 'CCATA': 'KL', 'CCATC': 'KM', 'CCATG': 'K4', 'CCATT': 'KN',
               'CCCAA': 'KO', 'CCCAC': 'K5', 'CCCAG': 'KP', 'CCCAT': 'KQ', 'CCCCA': 'K6', 'CCCCC': 'KR',
               'CCCCG': 'KS', 'CCCCT': 'K7', 'CCCGA': 'KU', 'CCCGC': 'KV', 'CCCGG': 'K8', 'CCCGT': 'KW',
               'CCCTA': 'KX', 'CCCTC': 'K9', 'CCCTG': 'KY', 'CCCTT': 'KZ', 'CCGAA': '3B', 'CCGAC': '3D',
               'CCGAG': '30', 'CCGAT': '3E', 'CCGCA': '3F', 'CCGCC': '31', 'CCGCG': '3H', 'CCGCT': '3I',
               'CCGGA': '32', 'CCGGC': '3J', 'CCGGG': '3K', 'CCGGT': '33', 'CCGTA': '3L', 'CCGTC': '3M',
               'CCGTG': '34', 'CCGTT': '3N', 'CCTAA': '3O', 'CCTAC': '35', 'CCTAG': '3P', 'CCTAT': '3Q',
               'CCTCA': '36', 'CCTCC': '3R', 'CCTCG': '3S', 'CCTCT': '37', 'CCTGA': '3U', 'CCTGC': '3V',
               'CCTGG': '38', 'CCTGT': '3W', 'CCTTA': '3X', 'CCTTC': '39', 'CCTTG': '3Y', 'CCTTT': '3Z',
               'CGAAA': 'LB', 'CGAAC': 'LD', 'CGAAG': 'L0', 'CGAAT': 'LE', 'CGACA': 'LF', 'CGACC': 'L1',
               'CGACG': 'LH', 'CGACT': 'LI', 'CGAGA': 'L2', 'CGAGC': 'LJ', 'CGAGG': 'LK', 'CGAGT': 'L3',
               'CGATA': 'LL', 'CGATC': 'LM', 'CGATG': 'L4', 'CGATT': 'LN', 'CGCAA': 'LO', 'CGCAC': 'L5',
               'CGCAG': 'LP', 'CGCAT': 'LQ', 'CGCCA': 'L6', 'CGCCC': 'LR', 'CGCCG': 'LS', 'CGCCT': 'L7',
               'CGCGA': 'LU', 'CGCGC': 'LV', 'CGCGG': 'L8', 'CGCGT': 'LW', 'CGCTA': 'LX', 'CGCTC': 'L9',
               'CGCTG': 'LY', 'CGCTT': 'LZ', 'CGGAA': 'MB', 'CGGAC': 'MD', 'CGGAG': 'M0', 'CGGAT': 'ME',
               'CGGCA': 'MF', 'CGGCC': 'M1', 'CGGCG': 'MH', 'CGGCT': 'MI', 'CGGGA': 'M2', 'CGGGC': 'MJ',
               'CGGGG': 'MK', 'CGGGT': 'M3', 'CGGTA': 'ML', 'CGGTC': 'MM', 'CGGTG': 'M4', 'CGGTT': 'MN',
               'CGTAA': 'MO', 'CGTAC': 'M5', 'CGTAG': 'MP', 'CGTAT': 'MQ', 'CGTCA': 'M6', 'CGTCC': 'MR',
               'CGTCG': 'MS', 'CGTCT': 'M7', 'CGTGA': 'MU', 'CGTGC': 'MV', 'CGTGG': 'M8', 'CGTGT': 'MW',
               'CGTTA': 'MX', 'CGTTC': 'M9', 'CGTTG': 'MY', 'CGTTT': 'MZ', 'CTAAA': '4B', 'CTAAC': '4D',
               'CTAAG': '40', 'CTAAT': '4E', 'CTACA': '4F', 'CTACC': '41', 'CTACG': '4H', 'CTACT': '4I',
               'CTAGA': '42', 'CTAGC': '4J', 'CTAGG': '4K', 'CTAGT': '43', 'CTATA': '4L', 'CTATC': '4M',
               'CTATG': '44', 'CTATT': '4N', 'CTCAA': '4O', 'CTCAC': '45', 'CTCAG': '4P', 'CTCAT': '4Q',
               'CTCCA': '46', 'CTCCC': '4R', 'CTCCG': '4S', 'CTCCT': '47', 'CTCGA': '4U', 'CTCGC': '4V',
               'CTCGG': '48', 'CTCGT': '4W', 'CTCTA': '4X', 'CTCTC': '49', 'CTCTG': '4Y', 'CTCTT': '4Z',
               'CTGAA': 'NB', 'CTGAC': 'ND', 'CTGAG': 'N0', 'CTGAT': 'NE', 'CTGCA': 'NF', 'CTGCC': 'N1',
               'CTGCG': 'NH', 'CTGCT': 'NI', 'CTGGA': 'N2', 'CTGGC': 'NJ', 'CTGGG': 'NK', 'CTGGT': 'N3',
               'CTGTA': 'NL', 'CTGTC': 'NM', 'CTGTG': 'N4', 'CTGTT': 'NN', 'CTTAA': 'NO', 'CTTAC': 'N5',
               'CTTAG': 'NP', 'CTTAT': 'NQ', 'CTTCA': 'N6', 'CTTCC': 'NR', 'CTTCG': 'NS', 'CTTCT': 'N7',
               'CTTGA': 'NU', 'CTTGC': 'NV', 'CTTGG': 'N8', 'CTTGT': 'NW', 'CTTTA': 'NX', 'CTTTC': 'N9',
               'CTTTG': 'NY', 'CTTTT': 'NZ', 'GAAAA': 'OB', 'GAAAC': 'OD', 'GAAAG': 'O0', 'GAAAT': 'OE',
               'GAACA': 'OF', 'GAACC': 'O1', 'GAACG': 'OH', 'GAACT': 'OI', 'GAAGA': 'O2', 'GAAGC': 'OJ',
               'GAAGG': 'OK', 'GAAGT': 'O3', 'GAATA': 'OL', 'GAATC': 'OM', 'GAATG': 'O4', 'GAATT': 'ON',
               'GACAA': 'OO', 'GACAC': 'O5', 'GACAG': 'OP', 'GACAT': 'OQ', 'GACCA': 'O6', 'GACCC': 'OR',
               'GACCG': 'OS', 'GACCT': 'O7', 'GACGA': 'OU', 'GACGC': 'OV', 'GACGG': 'O8', 'GACGT': 'OW',
               'GACTA': 'OX', 'GACTC': 'O9', 'GACTG': 'OY', 'GACTT': 'OZ', 'GAGAA': '5B', 'GAGAC': '5D',
               'GAGAG': '50', 'GAGAT': '5E', 'GAGCA': '5F', 'GAGCC': '51', 'GAGCG': '5H', 'GAGCT': '5I',
               'GAGGA': '52', 'GAGGC': '5J', 'GAGGG': '5K', 'GAGGT': '53', 'GAGTA': '5L', 'GAGTC': '5M',
               'GAGTG': '54', 'GAGTT': '5N', 'GATAA': '5O', 'GATAC': '55', 'GATAG': '5P', 'GATAT': '5Q',
               'GATCA': '56', 'GATCC': '5R', 'GATCG': '5S', 'GATCT': '57', 'GATGA': '5U', 'GATGC': '5V',
               'GATGG': '58', 'GATGT': '5W', 'GATTA': '5X', 'GATTC': '59', 'GATTG': '5Y', 'GATTT': '5Z',
               'GCAAA': 'PB', 'GCAAC': 'PD', 'GCAAG': 'P0', 'GCAAT': 'PE', 'GCACA': 'PF', 'GCACC': 'P1',
               'GCACG': 'PH', 'GCACT': 'PI', 'GCAGA': 'P2', 'GCAGC': 'PJ', 'GCAGG': 'PK', 'GCAGT': 'P3',
               'GCATA': 'PL', 'GCATC': 'PM', 'GCATG': 'P4', 'GCATT': 'PN', 'GCCAA': 'PO', 'GCCAC': 'P5',
               'GCCAG': 'PP', 'GCCAT': 'PQ', 'GCCCA': 'P6', 'GCCCC': 'PR', 'GCCCG': 'PS', 'GCCCT': 'P7',
               'GCCGA': 'PU', 'GCCGC': 'PV', 'GCCGG': 'P8', 'GCCGT': 'PW', 'GCCTA': 'PX', 'GCCTC': 'P9',
               'GCCTG': 'PY', 'GCCTT': 'PZ', 'GCGAA': 'QB', 'GCGAC': 'QD', 'GCGAG': 'Q0', 'GCGAT': 'QE',
               'GCGCA': 'QF', 'GCGCC': 'Q1', 'GCGCG': 'QH', 'GCGCT': 'QI', 'GCGGA': 'Q2', 'GCGGC': 'QJ',
               'GCGGG': 'QK', 'GCGGT': 'Q3', 'GCGTA': 'QL', 'GCGTC': 'QM', 'GCGTG': 'Q4', 'GCGTT': 'QN',
               'GCTAA': 'QO', 'GCTAC': 'Q5', 'GCTAG': 'QP', 'GCTAT': 'QQ', 'GCTCA': 'Q6', 'GCTCC': 'QR',
               'GCTCG': 'QS', 'GCTCT': 'Q7', 'GCTGA': 'QU', 'GCTGC': 'QV', 'GCTGG': 'Q8', 'GCTGT': 'QW',
               'GCTTA': 'QX', 'GCTTC': 'Q9', 'GCTTG': 'QY', 'GCTTT': 'QZ', 'GGAAA': '6B', 'GGAAC': '6D',
               'GGAAG': '60', 'GGAAT': '6E', 'GGACA': '6F', 'GGACC': '61', 'GGACG': '6H', 'GGACT': '6I',
               'GGAGA': '62', 'GGAGC': '6J', 'GGAGG': '6K', 'GGAGT': '63', 'GGATA': '6L', 'GGATC': '6M',
               'GGATG': '64', 'GGATT': '6N', 'GGCAA': '6O', 'GGCAC': '65', 'GGCAG': '6P', 'GGCAT': '6Q',
               'GGCCA': '66', 'GGCCC': '6R', 'GGCCG': '6S', 'GGCCT': '67', 'GGCGA': '6U', 'GGCGC': '6V',
               'GGCGG': '68', 'GGCGT': '6W', 'GGCTA': '6X', 'GGCTC': '69', 'GGCTG': '6Y', 'GGCTT': '6Z',
               'GGGAA': 'RB', 'GGGAC': 'RD', 'GGGAG': 'R0', 'GGGAT': 'RE', 'GGGCA': 'RF', 'GGGCC': 'R1',
               'GGGCG': 'RH', 'GGGCT': 'RI', 'GGGGA': 'R2', 'GGGGC': 'RJ', 'GGGGG': 'RK', 'GGGGT': 'R3',
               'GGGTA': 'RL', 'GGGTC': 'RM', 'GGGTG': 'R4', 'GGGTT': 'RN', 'GGTAA': 'RO', 'GGTAC': 'R5',
               'GGTAG': 'RP', 'GGTAT': 'RQ', 'GGTCA': 'R6', 'GGTCC': 'RR', 'GGTCG': 'RS', 'GGTCT': 'R7',
               'GGTGA': 'RU', 'GGTGC': 'RV', 'GGTGG': 'R8', 'GGTGT': 'RW', 'GGTTA': 'RX', 'GGTTC': 'R9',
               'GGTTG': 'RY', 'GGTTT': 'RZ', 'GTAAA': 'SB', 'GTAAC': 'SD', 'GTAAG': 'S0', 'GTAAT': 'SE',
               'GTACA': 'SF', 'GTACC': 'S1', 'GTACG': 'SH', 'GTACT': 'SI', 'GTAGA': 'S2', 'GTAGC': 'SJ',
               'GTAGG': 'SK', 'GTAGT': 'S3', 'GTATA': 'SL', 'GTATC': 'SM', 'GTATG': 'S4', 'GTATT': 'SN',
               'GTCAA': 'SO', 'GTCAC': 'S5', 'GTCAG': 'SP', 'GTCAT': 'SQ', 'GTCCA': 'S6', 'GTCCC': 'SR',
               'GTCCG': 'SS', 'GTCCT': 'S7', 'GTCGA': 'SU', 'GTCGC': 'SV', 'GTCGG': 'S8', 'GTCGT': 'SW',
               'GTCTA': 'SX', 'GTCTC': 'S9', 'GTCTG': 'SY', 'GTCTT': 'SZ', 'GTGAA': '7B', 'GTGAC': '7D',
               'GTGAG': '70', 'GTGAT': '7E', 'GTGCA': '7F', 'GTGCC': '71', 'GTGCG': '7H', 'GTGCT': '7I',
               'GTGGA': '72', 'GTGGC': '7J', 'GTGGG': '7K', 'GTGGT': '73', 'GTGTA': '7L', 'GTGTC': '7M',
               'GTGTG': '74', 'GTGTT': '7N', 'GTTAA': '7O', 'GTTAC': '75', 'GTTAG': '7P', 'GTTAT': '7Q',
               'GTTCA': '76', 'GTTCC': '7R', 'GTTCG': '7S', 'GTTCT': '77', 'GTTGA': '7U', 'GTTGC': '7V',
               'GTTGG': '78', 'GTTGT': '7W', 'GTTTA': '7X', 'GTTTC': '79', 'GTTTG': '7Y', 'GTTTT': '7Z',
               'TAAAA': 'UB', 'TAAAC': 'UD', 'TAAAG': 'U0', 'TAAAT': 'UE', 'TAACA': 'UF', 'TAACC': 'U1',
               'TAACG': 'UH', 'TAACT': 'UI', 'TAAGA': 'U2', 'TAAGC': 'UJ', 'TAAGG': 'UK', 'TAAGT': 'U3',
               'TAATA': 'UL', 'TAATC': 'UM', 'TAATG': 'U4', 'TAATT': 'UN', 'TACAA': 'UO', 'TACAC': 'U5',
               'TACAG': 'UP', 'TACAT': 'UQ', 'TACCA': 'U6', 'TACCC': 'UR', 'TACCG': 'US', 'TACCT': 'U7',
               'TACGA': 'UU', 'TACGC': 'UV', 'TACGG': 'U8', 'TACGT': 'UW', 'TACTA': 'UX', 'TACTC': 'U9',
               'TACTG': 'UY', 'TACTT': 'UZ', 'TAGAA': 'VB', 'TAGAC': 'VD', 'TAGAG': 'V0', 'TAGAT': 'VE',
               'TAGCA': 'VF', 'TAGCC': 'V1', 'TAGCG': 'VH', 'TAGCT': 'VI', 'TAGGA': 'V2', 'TAGGC': 'VJ',
               'TAGGG': 'VK', 'TAGGT': 'V3', 'TAGTA': 'VL', 'TAGTC': 'VM', 'TAGTG': 'V4', 'TAGTT': 'VN',
               'TATAA': 'VO', 'TATAC': 'V5', 'TATAG': 'VP', 'TATAT': 'VQ', 'TATCA': 'V6', 'TATCC': 'VR',
               'TATCG': 'VS', 'TATCT': 'V7', 'TATGA': 'VU', 'TATGC': 'VV', 'TATGG': 'V8', 'TATGT': 'VW',
               'TATTA': 'VX', 'TATTC': 'V9', 'TATTG': 'VY', 'TATTT': 'VZ', 'TCAAA': '8B', 'TCAAC': '8D',
               'TCAAG': '80', 'TCAAT': '8E', 'TCACA': '8F', 'TCACC': '81', 'TCACG': '8H', 'TCACT': '8I',
               'TCAGA': '82', 'TCAGC': '8J', 'TCAGG': '8K', 'TCAGT': '83', 'TCATA': '8L', 'TCATC': '8M',
               'TCATG': '84', 'TCATT': '8N', 'TCCAA': '8O', 'TCCAC': '85', 'TCCAG': '8P', 'TCCAT': '8Q',
               'TCCCA': '86', 'TCCCC': '8R', 'TCCCG': '8S', 'TCCCT': '87', 'TCCGA': '8U', 'TCCGC': '8V',
               'TCCGG': '88', 'TCCGT': '8W', 'TCCTA': '8X', 'TCCTC': '89', 'TCCTG': '8Y', 'TCCTT': '8Z',
               'TCGAA': 'WB', 'TCGAC': 'WD', 'TCGAG': 'W0', 'TCGAT': 'WE', 'TCGCA': 'WF', 'TCGCC': 'W1',
               'TCGCG': 'WH', 'TCGCT': 'WI', 'TCGGA': 'W2', 'TCGGC': 'WJ', 'TCGGG': 'WK', 'TCGGT': 'W3',
               'TCGTA': 'WL', 'TCGTC': 'WM', 'TCGTG': 'W4', 'TCGTT': 'WN', 'TCTAA': 'WO', 'TCTAC': 'W5',
               'TCTAG': 'WP', 'TCTAT': 'WQ', 'TCTCA': 'W6', 'TCTCC': 'WR', 'TCTCG': 'WS', 'TCTCT': 'W7',
               'TCTGA': 'WU', 'TCTGC': 'WV', 'TCTGG': 'W8', 'TCTGT': 'WW', 'TCTTA': 'WX', 'TCTTC': 'W9',
               'TCTTG': 'WY', 'TCTTT': 'WZ', 'TGAAA': 'XB', 'TGAAC': 'XD', 'TGAAG': 'X0', 'TGAAT': 'XE',
               'TGACA': 'XF', 'TGACC': 'X1', 'TGACG': 'XH', 'TGACT': 'XI', 'TGAGA': 'X2', 'TGAGC': 'XJ',
               'TGAGG': 'XK', 'TGAGT': 'X3', 'TGATA': 'XL', 'TGATC': 'XM', 'TGATG': 'X4', 'TGATT': 'XN',
               'TGCAA': 'XO', 'TGCAC': 'X5', 'TGCAG': 'XP', 'TGCAT': 'XQ', 'TGCCA': 'X6', 'TGCCC': 'XR',
               'TGCCG': 'XS', 'TGCCT': 'X7', 'TGCGA': 'XU', 'TGCGC': 'XV', 'TGCGG': 'X8', 'TGCGT': 'XW',
               'TGCTA': 'XX', 'TGCTC': 'X9', 'TGCTG': 'XY', 'TGCTT': 'XZ', 'TGGAA': '9B', 'TGGAC': '9D',
               'TGGAG': '90', 'TGGAT': '9E', 'TGGCA': '9F', 'TGGCC': '91', 'TGGCG': '9H', 'TGGCT': '9I',
               'TGGGA': '92', 'TGGGC': '9J', 'TGGGG': '9K', 'TGGGT': '93', 'TGGTA': '9L', 'TGGTC': '9M',
               'TGGTG': '94', 'TGGTT': '9N', 'TGTAA': '9O', 'TGTAC': '95', 'TGTAG': '9P', 'TGTAT': '9Q',
               'TGTCA': '96', 'TGTCC': '9R', 'TGTCG': '9S', 'TGTCT': '97', 'TGTGA': '9U', 'TGTGC': '9V',
               'TGTGG': '98', 'TGTGT': '9W', 'TGTTA': '9X', 'TGTTC': '99', 'TGTTG': '9Y', 'TGTTT': '9Z',
               'TTAAA': 'YB', 'TTAAC': 'YD', 'TTAAG': 'Y0', 'TTAAT': 'YE', 'TTACA': 'YF', 'TTACC': 'Y1',
               'TTACG': 'YH', 'TTACT': 'YI', 'TTAGA': 'Y2', 'TTAGC': 'YJ', 'TTAGG': 'YK', 'TTAGT': 'Y3',
               'TTATA': 'YL', 'TTATC': 'YM', 'TTATG': 'Y4', 'TTATT': 'YN', 'TTCAA': 'YO', 'TTCAC': 'Y5',
               'TTCAG': 'YP', 'TTCAT': 'YQ', 'TTCCA': 'Y6', 'TTCCC': 'YR', 'TTCCG': 'YS', 'TTCCT': 'Y7',
               'TTCGA': 'YU', 'TTCGC': 'YV', 'TTCGG': 'Y8', 'TTCGT': 'YW', 'TTCTA': 'YX', 'TTCTC': 'Y9',
               'TTCTG': 'YY', 'TTCTT': 'YZ', 'TTGAA': 'ZB', 'TTGAC': 'ZD', 'TTGAG': 'Z0', 'TTGAT': 'ZE',
               'TTGCA': 'ZF', 'TTGCC': 'Z1', 'TTGCG': 'ZH', 'TTGCT': 'ZI', 'TTGGA': 'Z2', 'TTGGC': 'ZJ',
               'TTGGG': 'ZK', 'TTGGT': 'Z3', 'TTGTA': 'ZL', 'TTGTC': 'ZM', 'TTGTG': 'Z4', 'TTGTT': 'ZN',
               'TTTAA': 'ZO', 'TTTAC': 'Z5', 'TTTAG': 'ZP', 'TTTAT': 'ZQ', 'TTTCA': 'Z6', 'TTTCC': 'ZR',
               'TTTCG': 'ZS', 'TTTCT': 'Z7', 'TTTGA': 'ZU', 'TTTGC': 'ZV', 'TTTGG': 'Z8', 'TTTGT': 'ZW',
               'TTTTA': 'ZX', 'TTTTC': 'Z9', 'TTTTG': 'ZY', 'TTTTT': 'ZZ', 'A': 'B', 'C': 'D',
               'G': '0', 'T': 'E', 'AA': 'F', 'AC': '1', 'AG': 'H', 'AT': 'I',
               'CA': '2', 'CC': 'J', 'CG': 'K', 'CT': '3', 'GA': 'L', 'GC': 'M',
               'GG': '4', 'GT': 'N', 'TA': 'O', 'TC': '5', 'TG': 'P', 'TT': 'Q',
               'AAA': '6', 'AAC': 'R', 'AAG': 'S', 'AAT': '7', 'ACA': 'U', 'ACC': 'V',
               'ACG': '8', 'ACT': 'W', 'AGA': 'X', 'AGC': '9', 'AGG': 'Y', 'AGT': 'Z',
               'ATA': 'DB', 'ATC': 'DD', 'ATG': 'D0', 'ATT': 'DE', 'CAA': 'DF', 'CAC': 'D1',
               'CAG': 'DH', 'CAT': 'DI', 'CCA': 'D2', 'CCC': 'DJ', 'CCG': 'DK', 'CCT': 'D3',
               'CGA': 'DL', 'CGC': 'DM', 'CGG': 'D4', 'CGT': 'DN', 'CTA': 'DO', 'CTC': 'D5',
               'CTG': 'DP', 'CTT': 'DQ', 'GAA': 'D6', 'GAC': 'DR', 'GAG': 'DS', 'GAT': 'D7',
               'GCA': 'DU', 'GCC': 'DV', 'GCG': 'D8', 'GCT': 'DW', 'GGA': 'DX', 'GGC': 'D9',
               'GGG': 'DY', 'GGT': 'DZ', 'GTA': '0B', 'GTC': '0D', 'GTG': '00', 'GTT': '0E',
               'TAA': '0F', 'TAC': '01', 'TAG': '0H', 'TAT': '0I', 'TCA': '02', 'TCC': '0J',
               'TCG': '0K', 'TCT': '03', 'TGA': '0L', 'TGC': '0M', 'TGG': '04', 'TGT': '0N',
               'TTA': '0O', 'TTC': '05', 'TTG': '0P', 'TTT': '0Q', 'AAAA': '06', 'AAAC': '0R',
               'AAAG': '0S', 'AAAT': '07', 'AACA': '0U', 'AACC': '0V', 'AACG': '08', 'AACT': '0W',
               'AAGA': '0X', 'AAGC': '09', 'AAGG': '0Y', 'AAGT': '0Z', 'AATA': 'EB', 'AATC': 'ED',
               'AATG': 'E0', 'AATT': 'EE', 'ACAA': 'EF', 'ACAC': 'E1', 'ACAG': 'EH', 'ACAT': 'EI',
               'ACCA': 'E2', 'ACCC': 'EJ', 'ACCG': 'EK', 'ACCT': 'E3', 'ACGA': 'EL', 'ACGC': 'EM',
               'ACGG': 'E4', 'ACGT': 'EN', 'ACTA': 'EO', 'ACTC': 'E5', 'ACTG': 'EP', 'ACTT': 'EQ',
               'AGAA': 'E6', 'AGAC': 'ER', 'AGAG': 'ES', 'AGAT': 'E7', 'AGCA': 'EU', 'AGCC': 'EV',
               'AGCG': 'E8', 'AGCT': 'EW', 'AGGA': 'EX', 'AGGC': 'E9', 'AGGG': 'EY', 'AGGT': 'EZ',
               'AGTA': 'FB', 'AGTC': 'FD', 'AGTG': 'F0', 'AGTT': 'FE', 'ATAA': 'FF', 'ATAC': 'F1',
               'ATAG': 'FH', 'ATAT': 'FI', 'ATCA': 'F2', 'ATCC': 'FJ', 'ATCG': 'FK', 'ATCT': 'F3',
               'ATGA': 'FL', 'ATGC': 'FM', 'ATGG': 'F4', 'ATGT': 'FN', 'ATTA': 'FO', 'ATTC': 'F5',
               'ATTG': 'FP', 'ATTT': 'FQ', 'CAAA': 'F6', 'CAAC': 'FR', 'CAAG': 'FS', 'CAAT': 'F7',
               'CACA': 'FU', 'CACC': 'FV', 'CACG': 'F8', 'CACT': 'FW', 'CAGA': 'FX', 'CAGC': 'F9',
               'CAGG': 'FY', 'CAGT': 'FZ', 'CATA': '1B', 'CATC': '1D', 'CATG': '10', 'CATT': '1E',
               'CCAA': '1F', 'CCAC': '11', 'CCAG': '1H', 'CCAT': '1I', 'CCCA': '12', 'CCCC': '1J',
               'CCCG': '1K', 'CCCT': '13', 'CCGA': '1L', 'CCGC': '1M', 'CCGG': '14', 'CCGT': '1N',
               'CCTA': '1O', 'CCTC': '15', 'CCTG': '1P', 'CCTT': '1Q', 'CGAA': '16', 'CGAC': '1R',
               'CGAG': '1S', 'CGAT': '17', 'CGCA': '1U', 'CGCC': '1V', 'CGCG': '18', 'CGCT': '1W',
               'CGGA': '1X', 'CGGC': '19', 'CGGG': '1Y', 'CGGT': '1Z', 'CGTA': 'HB', 'CGTC': 'HD',
               'CGTG': 'H0', 'CGTT': 'HE', 'CTAA': 'HF', 'CTAC': 'H1', 'CTAG': 'HH', 'CTAT': 'HI',
               'CTCA': 'H2', 'CTCC': 'HJ', 'CTCG': 'HK', 'CTCT': 'H3', 'CTGA': 'HL', 'CTGC': 'HM',
               'CTGG': 'H4', 'CTGT': 'HN', 'CTTA': 'HO', 'CTTC': 'H5', 'CTTG': 'HP', 'CTTT': 'HQ',
               'GAAA': 'H6', 'GAAC': 'HR', 'GAAG': 'HS', 'GAAT': 'H7', 'GACA': 'HU', 'GACC': 'HV',
               'GACG': 'H8', 'GACT': 'HW', 'GAGA': 'HX', 'GAGC': 'H9', 'GAGG': 'HY', 'GAGT': 'HZ',
               'GATA': 'IB', 'GATC': 'ID', 'GATG': 'I0', 'GATT': 'IE', 'GCAA': 'IF', 'GCAC': 'I1',
               'GCAG': 'IH', 'GCAT': 'II', 'GCCA': 'I2', 'GCCC': 'IJ', 'GCCG': 'IK', 'GCCT': 'I3',
               'GCGA': 'IL', 'GCGC': 'IM', 'GCGG': 'I4', 'GCGT': 'IN', 'GCTA': 'IO', 'GCTC': 'I5',
               'GCTG': 'IP', 'GCTT': 'IQ', 'GGAA': 'I6', 'GGAC': 'IR', 'GGAG': 'IS', 'GGAT': 'I7',
               'GGCA': 'IU', 'GGCC': 'IV', 'GGCG': 'I8', 'GGCT': 'IW', 'GGGA': 'IX', 'GGGC': 'I9',
               'GGGG': 'IY', 'GGGT': 'IZ', 'GGTA': '2B', 'GGTC': '2D', 'GGTG': '20', 'GGTT': '2E',
               'GTAA': '2F', 'GTAC': '21', 'GTAG': '2H', 'GTAT': '2I', 'GTCA': '22', 'GTCC': '2J',
               'GTCG': '2K', 'GTCT': '23', 'GTGA': '2L', 'GTGC': '2M', 'GTGG': '24', 'GTGT': '2N',
               'GTTA': '2O', 'GTTC': '25', 'GTTG': '2P', 'GTTT': '2Q', 'TAAA': '26', 'TAAC': '2R',
               'TAAG': '2S', 'TAAT': '27', 'TACA': '2U', 'TACC': '2V', 'TACG': '28', 'TACT': '2W',
               'TAGA': '2X', 'TAGC': '29', 'TAGG': '2Y', 'TAGT': '2Z', 'TATA': 'JB', 'TATC': 'JD',
               'TATG': 'J0', 'TATT': 'JE', 'TCAA': 'JF', 'TCAC': 'J1', 'TCAG': 'JH', 'TCAT': 'JI',
               'TCCA': 'J2', 'TCCC': 'JJ', 'TCCG': 'JK', 'TCCT': 'J3', 'TCGA': 'JL', 'TCGC': 'JM',
               'TCGG': 'J4', 'TCGT': 'JN', 'TCTA': 'JO', 'TCTC': 'J5', 'TCTG': 'JP', 'TCTT': 'JQ',
               'TGAA': 'J6', 'TGAC': 'JR', 'TGAG': 'JS', 'TGAT': 'J7', 'TGCA': 'JU', 'TGCC': 'JV',
               'TGCG': 'J8', 'TGCT': 'JW', 'TGGA': 'JX', 'TGGC': 'J9', 'TGGG': 'JY', 'TGGT': 'JZ',
               'TGTA': 'KB', 'TGTC': 'KD', 'TGTG': 'K0', 'TGTT': 'KE', 'TTAA': 'KF', 'TTAC': 'K1',
               'TTAG': 'KH', 'TTAT': 'KI', 'TTCA': 'K2', 'TTCC': 'KJ', 'TTCG': 'KK', 'TTCT': 'K3',
               'TTGA': 'KL', 'TTGC': 'KM', 'TTGG': 'K4', 'TTGT': 'KN', 'TTTA': 'KO', 'TTTC': 'K5',
               'TTTG': 'KP', 'TTTT': 'KQ'}

decode_hash = {'BB-5': 'AAAAA', 'BD-5': 'AAAAC', 'B0-5': 'AAAAG', 'BE-5': 'AAAAT', 'BF-5': 'AAACA', 'B1-5': 'AAACC',
               'BH-5': 'AAACG', 'BI-5': 'AAACT', 'B2-5': 'AAAGA', 'BJ-5': 'AAAGC', 'BK-5': 'AAAGG', 'B3-5': 'AAAGT',
               'BL-5': 'AAATA', 'BM-5': 'AAATC', 'B4-5': 'AAATG', 'BN-5': 'AAATT', 'BO-5': 'AACAA', 'B5-5': 'AACAC',
               'BP-5': 'AACAG', 'BQ-5': 'AACAT', 'B6-5': 'AACCA', 'BR-5': 'AACCC', 'BS-5': 'AACCG', 'B7-5': 'AACCT',
               'BU-5': 'AACGA', 'BV-5': 'AACGC', 'B8-5': 'AACGG', 'BW-5': 'AACGT', 'BX-5': 'AACTA', 'B9-5': 'AACTC',
               'BY-5': 'AACTG', 'BZ-5': 'AACTT', 'DB-5': 'AAGAA', 'DD-5': 'AAGAC', 'D0-5': 'AAGAG', 'DE-5': 'AAGAT',
               'DF-5': 'AAGCA', 'D1-5': 'AAGCC', 'DH-5': 'AAGCG', 'DI-5': 'AAGCT', 'D2-5': 'AAGGA', 'DJ-5': 'AAGGC',
               'DK-5': 'AAGGG', 'D3-5': 'AAGGT', 'DL-5': 'AAGTA', 'DM-5': 'AAGTC', 'D4-5': 'AAGTG', 'DN-5': 'AAGTT',
               'DO-5': 'AATAA', 'D5-5': 'AATAC', 'DP-5': 'AATAG', 'DQ-5': 'AATAT', 'D6-5': 'AATCA', 'DR-5': 'AATCC',
               'DS-5': 'AATCG', 'D7-5': 'AATCT', 'DU-5': 'AATGA', 'DV-5': 'AATGC', 'D8-5': 'AATGG', 'DW-5': 'AATGT',
               'DX-5': 'AATTA', 'D9-5': 'AATTC', 'DY-5': 'AATTG', 'DZ-5': 'AATTT', '0B-5': 'ACAAA', '0D-5': 'ACAAC',
               '00-5': 'ACAAG', '0E-5': 'ACAAT', '0F-5': 'ACACA', '01-5': 'ACACC', '0H-5': 'ACACG', '0I-5': 'ACACT',
               '02-5': 'ACAGA', '0J-5': 'ACAGC', '0K-5': 'ACAGG', '03-5': 'ACAGT', '0L-5': 'ACATA', '0M-5': 'ACATC',
               '04-5': 'ACATG', '0N-5': 'ACATT', '0O-5': 'ACCAA', '05-5': 'ACCAC', '0P-5': 'ACCAG', '0Q-5': 'ACCAT',
               '06-5': 'ACCCA', '0R-5': 'ACCCC', '0S-5': 'ACCCG', '07-5': 'ACCCT', '0U-5': 'ACCGA', '0V-5': 'ACCGC',
               '08-5': 'ACCGG', '0W-5': 'ACCGT', '0X-5': 'ACCTA', '09-5': 'ACCTC', '0Y-5': 'ACCTG', '0Z-5': 'ACCTT',
               'EB-5': 'ACGAA', 'ED-5': 'ACGAC', 'E0-5': 'ACGAG', 'EE-5': 'ACGAT', 'EF-5': 'ACGCA', 'E1-5': 'ACGCC',
               'EH-5': 'ACGCG', 'EI-5': 'ACGCT', 'E2-5': 'ACGGA', 'EJ-5': 'ACGGC', 'EK-5': 'ACGGG', 'E3-5': 'ACGGT',
               'EL-5': 'ACGTA', 'EM-5': 'ACGTC', 'E4-5': 'ACGTG', 'EN-5': 'ACGTT', 'EO-5': 'ACTAA', 'E5-5': 'ACTAC',
               'EP-5': 'ACTAG', 'EQ-5': 'ACTAT', 'E6-5': 'ACTCA', 'ER-5': 'ACTCC', 'ES-5': 'ACTCG', 'E7-5': 'ACTCT',
               'EU-5': 'ACTGA', 'EV-5': 'ACTGC', 'E8-5': 'ACTGG', 'EW-5': 'ACTGT', 'EX-5': 'ACTTA', 'E9-5': 'ACTTC',
               'EY-5': 'ACTTG', 'EZ-5': 'ACTTT', 'FB-5': 'AGAAA', 'FD-5': 'AGAAC', 'F0-5': 'AGAAG', 'FE-5': 'AGAAT',
               'FF-5': 'AGACA', 'F1-5': 'AGACC', 'FH-5': 'AGACG', 'FI-5': 'AGACT', 'F2-5': 'AGAGA', 'FJ-5': 'AGAGC',
               'FK-5': 'AGAGG', 'F3-5': 'AGAGT', 'FL-5': 'AGATA', 'FM-5': 'AGATC', 'F4-5': 'AGATG', 'FN-5': 'AGATT',
               'FO-5': 'AGCAA', 'F5-5': 'AGCAC', 'FP-5': 'AGCAG', 'FQ-5': 'AGCAT', 'F6-5': 'AGCCA', 'FR-5': 'AGCCC',
               'FS-5': 'AGCCG', 'F7-5': 'AGCCT', 'FU-5': 'AGCGA', 'FV-5': 'AGCGC', 'F8-5': 'AGCGG', 'FW-5': 'AGCGT',
               'FX-5': 'AGCTA', 'F9-5': 'AGCTC', 'FY-5': 'AGCTG', 'FZ-5': 'AGCTT', '1B-5': 'AGGAA', '1D-5': 'AGGAC',
               '10-5': 'AGGAG', '1E-5': 'AGGAT', '1F-5': 'AGGCA', '11-5': 'AGGCC', '1H-5': 'AGGCG', '1I-5': 'AGGCT',
               '12-5': 'AGGGA', '1J-5': 'AGGGC', '1K-5': 'AGGGG', '13-5': 'AGGGT', '1L-5': 'AGGTA', '1M-5': 'AGGTC',
               '14-5': 'AGGTG', '1N-5': 'AGGTT', '1O-5': 'AGTAA', '15-5': 'AGTAC', '1P-5': 'AGTAG', '1Q-5': 'AGTAT',
               '16-5': 'AGTCA', '1R-5': 'AGTCC', '1S-5': 'AGTCG', '17-5': 'AGTCT', '1U-5': 'AGTGA', '1V-5': 'AGTGC',
               '18-5': 'AGTGG', '1W-5': 'AGTGT', '1X-5': 'AGTTA', '19-5': 'AGTTC', '1Y-5': 'AGTTG', '1Z-5': 'AGTTT',
               'HB-5': 'ATAAA', 'HD-5': 'ATAAC', 'H0-5': 'ATAAG', 'HE-5': 'ATAAT', 'HF-5': 'ATACA', 'H1-5': 'ATACC',
               'HH-5': 'ATACG', 'HI-5': 'ATACT', 'H2-5': 'ATAGA', 'HJ-5': 'ATAGC', 'HK-5': 'ATAGG', 'H3-5': 'ATAGT',
               'HL-5': 'ATATA', 'HM-5': 'ATATC', 'H4-5': 'ATATG', 'HN-5': 'ATATT', 'HO-5': 'ATCAA', 'H5-5': 'ATCAC',
               'HP-5': 'ATCAG', 'HQ-5': 'ATCAT', 'H6-5': 'ATCCA', 'HR-5': 'ATCCC', 'HS-5': 'ATCCG', 'H7-5': 'ATCCT',
               'HU-5': 'ATCGA', 'HV-5': 'ATCGC', 'H8-5': 'ATCGG', 'HW-5': 'ATCGT', 'HX-5': 'ATCTA', 'H9-5': 'ATCTC',
               'HY-5': 'ATCTG', 'HZ-5': 'ATCTT', 'IB-5': 'ATGAA', 'ID-5': 'ATGAC', 'I0-5': 'ATGAG', 'IE-5': 'ATGAT',
               'IF-5': 'ATGCA', 'I1-5': 'ATGCC', 'IH-5': 'ATGCG', 'II-5': 'ATGCT', 'I2-5': 'ATGGA', 'IJ-5': 'ATGGC',
               'IK-5': 'ATGGG', 'I3-5': 'ATGGT', 'IL-5': 'ATGTA', 'IM-5': 'ATGTC', 'I4-5': 'ATGTG', 'IN-5': 'ATGTT',
               'IO-5': 'ATTAA', 'I5-5': 'ATTAC', 'IP-5': 'ATTAG', 'IQ-5': 'ATTAT', 'I6-5': 'ATTCA', 'IR-5': 'ATTCC',
               'IS-5': 'ATTCG', 'I7-5': 'ATTCT', 'IU-5': 'ATTGA', 'IV-5': 'ATTGC', 'I8-5': 'ATTGG', 'IW-5': 'ATTGT',
               'IX-5': 'ATTTA', 'I9-5': 'ATTTC', 'IY-5': 'ATTTG', 'IZ-5': 'ATTTT', '2B-5': 'CAAAA', '2D-5': 'CAAAC',
               '20-5': 'CAAAG', '2E-5': 'CAAAT', '2F-5': 'CAACA', '21-5': 'CAACC', '2H-5': 'CAACG', '2I-5': 'CAACT',
               '22-5': 'CAAGA', '2J-5': 'CAAGC', '2K-5': 'CAAGG', '23-5': 'CAAGT', '2L-5': 'CAATA', '2M-5': 'CAATC',
               '24-5': 'CAATG', '2N-5': 'CAATT', '2O-5': 'CACAA', '25-5': 'CACAC', '2P-5': 'CACAG', '2Q-5': 'CACAT',
               '26-5': 'CACCA', '2R-5': 'CACCC', '2S-5': 'CACCG', '27-5': 'CACCT', '2U-5': 'CACGA', '2V-5': 'CACGC',
               '28-5': 'CACGG', '2W-5': 'CACGT', '2X-5': 'CACTA', '29-5': 'CACTC', '2Y-5': 'CACTG', '2Z-5': 'CACTT',
               'JB-5': 'CAGAA', 'JD-5': 'CAGAC', 'J0-5': 'CAGAG', 'JE-5': 'CAGAT', 'JF-5': 'CAGCA', 'J1-5': 'CAGCC',
               'JH-5': 'CAGCG', 'JI-5': 'CAGCT', 'J2-5': 'CAGGA', 'JJ-5': 'CAGGC', 'JK-5': 'CAGGG', 'J3-5': 'CAGGT',
               'JL-5': 'CAGTA', 'JM-5': 'CAGTC', 'J4-5': 'CAGTG', 'JN-5': 'CAGTT', 'JO-5': 'CATAA', 'J5-5': 'CATAC',
               'JP-5': 'CATAG', 'JQ-5': 'CATAT', 'J6-5': 'CATCA', 'JR-5': 'CATCC', 'JS-5': 'CATCG', 'J7-5': 'CATCT',
               'JU-5': 'CATGA', 'JV-5': 'CATGC', 'J8-5': 'CATGG', 'JW-5': 'CATGT', 'JX-5': 'CATTA', 'J9-5': 'CATTC',
               'JY-5': 'CATTG', 'JZ-5': 'CATTT', 'KB-5': 'CCAAA', 'KD-5': 'CCAAC', 'K0-5': 'CCAAG', 'KE-5': 'CCAAT',
               'KF-5': 'CCACA', 'K1-5': 'CCACC', 'KH-5': 'CCACG', 'KI-5': 'CCACT', 'K2-5': 'CCAGA', 'KJ-5': 'CCAGC',
               'KK-5': 'CCAGG', 'K3-5': 'CCAGT', 'KL-5': 'CCATA', 'KM-5': 'CCATC', 'K4-5': 'CCATG', 'KN-5': 'CCATT',
               'KO-5': 'CCCAA', 'K5-5': 'CCCAC', 'KP-5': 'CCCAG', 'KQ-5': 'CCCAT', 'K6-5': 'CCCCA', 'KR-5': 'CCCCC',
               'KS-5': 'CCCCG', 'K7-5': 'CCCCT', 'KU-5': 'CCCGA', 'KV-5': 'CCCGC', 'K8-5': 'CCCGG', 'KW-5': 'CCCGT',
               'KX-5': 'CCCTA', 'K9-5': 'CCCTC', 'KY-5': 'CCCTG', 'KZ-5': 'CCCTT', '3B-5': 'CCGAA', '3D-5': 'CCGAC',
               '30-5': 'CCGAG', '3E-5': 'CCGAT', '3F-5': 'CCGCA', '31-5': 'CCGCC', '3H-5': 'CCGCG', '3I-5': 'CCGCT',
               '32-5': 'CCGGA', '3J-5': 'CCGGC', '3K-5': 'CCGGG', '33-5': 'CCGGT', '3L-5': 'CCGTA', '3M-5': 'CCGTC',
               '34-5': 'CCGTG', '3N-5': 'CCGTT', '3O-5': 'CCTAA', '35-5': 'CCTAC', '3P-5': 'CCTAG', '3Q-5': 'CCTAT',
               '36-5': 'CCTCA', '3R-5': 'CCTCC', '3S-5': 'CCTCG', '37-5': 'CCTCT', '3U-5': 'CCTGA', '3V-5': 'CCTGC',
               '38-5': 'CCTGG', '3W-5': 'CCTGT', '3X-5': 'CCTTA', '39-5': 'CCTTC', '3Y-5': 'CCTTG', '3Z-5': 'CCTTT',
               'LB-5': 'CGAAA', 'LD-5': 'CGAAC', 'L0-5': 'CGAAG', 'LE-5': 'CGAAT', 'LF-5': 'CGACA', 'L1-5': 'CGACC',
               'LH-5': 'CGACG', 'LI-5': 'CGACT', 'L2-5': 'CGAGA', 'LJ-5': 'CGAGC', 'LK-5': 'CGAGG', 'L3-5': 'CGAGT',
               'LL-5': 'CGATA', 'LM-5': 'CGATC', 'L4-5': 'CGATG', 'LN-5': 'CGATT', 'LO-5': 'CGCAA', 'L5-5': 'CGCAC',
               'LP-5': 'CGCAG', 'LQ-5': 'CGCAT', 'L6-5': 'CGCCA', 'LR-5': 'CGCCC', 'LS-5': 'CGCCG', 'L7-5': 'CGCCT',
               'LU-5': 'CGCGA', 'LV-5': 'CGCGC', 'L8-5': 'CGCGG', 'LW-5': 'CGCGT', 'LX-5': 'CGCTA', 'L9-5': 'CGCTC',
               'LY-5': 'CGCTG', 'LZ-5': 'CGCTT', 'MB-5': 'CGGAA', 'MD-5': 'CGGAC', 'M0-5': 'CGGAG', 'ME-5': 'CGGAT',
               'MF-5': 'CGGCA', 'M1-5': 'CGGCC', 'MH-5': 'CGGCG', 'MI-5': 'CGGCT', 'M2-5': 'CGGGA', 'MJ-5': 'CGGGC',
               'MK-5': 'CGGGG', 'M3-5': 'CGGGT', 'ML-5': 'CGGTA', 'MM-5': 'CGGTC', 'M4-5': 'CGGTG', 'MN-5': 'CGGTT',
               'MO-5': 'CGTAA', 'M5-5': 'CGTAC', 'MP-5': 'CGTAG', 'MQ-5': 'CGTAT', 'M6-5': 'CGTCA', 'MR-5': 'CGTCC',
               'MS-5': 'CGTCG', 'M7-5': 'CGTCT', 'MU-5': 'CGTGA', 'MV-5': 'CGTGC', 'M8-5': 'CGTGG', 'MW-5': 'CGTGT',
               'MX-5': 'CGTTA', 'M9-5': 'CGTTC', 'MY-5': 'CGTTG', 'MZ-5': 'CGTTT', '4B-5': 'CTAAA', '4D-5': 'CTAAC',
               '40-5': 'CTAAG', '4E-5': 'CTAAT', '4F-5': 'CTACA', '41-5': 'CTACC', '4H-5': 'CTACG', '4I-5': 'CTACT',
               '42-5': 'CTAGA', '4J-5': 'CTAGC', '4K-5': 'CTAGG', '43-5': 'CTAGT', '4L-5': 'CTATA', '4M-5': 'CTATC',
               '44-5': 'CTATG', '4N-5': 'CTATT', '4O-5': 'CTCAA', '45-5': 'CTCAC', '4P-5': 'CTCAG', '4Q-5': 'CTCAT',
               '46-5': 'CTCCA', '4R-5': 'CTCCC', '4S-5': 'CTCCG', '47-5': 'CTCCT', '4U-5': 'CTCGA', '4V-5': 'CTCGC',
               '48-5': 'CTCGG', '4W-5': 'CTCGT', '4X-5': 'CTCTA', '49-5': 'CTCTC', '4Y-5': 'CTCTG', '4Z-5': 'CTCTT',
               'NB-5': 'CTGAA', 'ND-5': 'CTGAC', 'N0-5': 'CTGAG', 'NE-5': 'CTGAT', 'NF-5': 'CTGCA', 'N1-5': 'CTGCC',
               'NH-5': 'CTGCG', 'NI-5': 'CTGCT', 'N2-5': 'CTGGA', 'NJ-5': 'CTGGC', 'NK-5': 'CTGGG', 'N3-5': 'CTGGT',
               'NL-5': 'CTGTA', 'NM-5': 'CTGTC', 'N4-5': 'CTGTG', 'NN-5': 'CTGTT', 'NO-5': 'CTTAA', 'N5-5': 'CTTAC',
               'NP-5': 'CTTAG', 'NQ-5': 'CTTAT', 'N6-5': 'CTTCA', 'NR-5': 'CTTCC', 'NS-5': 'CTTCG', 'N7-5': 'CTTCT',
               'NU-5': 'CTTGA', 'NV-5': 'CTTGC', 'N8-5': 'CTTGG', 'NW-5': 'CTTGT', 'NX-5': 'CTTTA', 'N9-5': 'CTTTC',
               'NY-5': 'CTTTG', 'NZ-5': 'CTTTT', 'OB-5': 'GAAAA', 'OD-5': 'GAAAC', 'O0-5': 'GAAAG', 'OE-5': 'GAAAT',
               'OF-5': 'GAACA', 'O1-5': 'GAACC', 'OH-5': 'GAACG', 'OI-5': 'GAACT', 'O2-5': 'GAAGA', 'OJ-5': 'GAAGC',
               'OK-5': 'GAAGG', 'O3-5': 'GAAGT', 'OL-5': 'GAATA', 'OM-5': 'GAATC', 'O4-5': 'GAATG', 'ON-5': 'GAATT',
               'OO-5': 'GACAA', 'O5-5': 'GACAC', 'OP-5': 'GACAG', 'OQ-5': 'GACAT', 'O6-5': 'GACCA', 'OR-5': 'GACCC',
               'OS-5': 'GACCG', 'O7-5': 'GACCT', 'OU-5': 'GACGA', 'OV-5': 'GACGC', 'O8-5': 'GACGG', 'OW-5': 'GACGT',
               'OX-5': 'GACTA', 'O9-5': 'GACTC', 'OY-5': 'GACTG', 'OZ-5': 'GACTT', '5B-5': 'GAGAA', '5D-5': 'GAGAC',
               '50-5': 'GAGAG', '5E-5': 'GAGAT', '5F-5': 'GAGCA', '51-5': 'GAGCC', '5H-5': 'GAGCG', '5I-5': 'GAGCT',
               '52-5': 'GAGGA', '5J-5': 'GAGGC', '5K-5': 'GAGGG', '53-5': 'GAGGT', '5L-5': 'GAGTA', '5M-5': 'GAGTC',
               '54-5': 'GAGTG', '5N-5': 'GAGTT', '5O-5': 'GATAA', '55-5': 'GATAC', '5P-5': 'GATAG', '5Q-5': 'GATAT',
               '56-5': 'GATCA', '5R-5': 'GATCC', '5S-5': 'GATCG', '57-5': 'GATCT', '5U-5': 'GATGA', '5V-5': 'GATGC',
               '58-5': 'GATGG', '5W-5': 'GATGT', '5X-5': 'GATTA', '59-5': 'GATTC', '5Y-5': 'GATTG', '5Z-5': 'GATTT',
               'PB-5': 'GCAAA', 'PD-5': 'GCAAC', 'P0-5': 'GCAAG', 'PE-5': 'GCAAT', 'PF-5': 'GCACA', 'P1-5': 'GCACC',
               'PH-5': 'GCACG', 'PI-5': 'GCACT', 'P2-5': 'GCAGA', 'PJ-5': 'GCAGC', 'PK-5': 'GCAGG', 'P3-5': 'GCAGT',
               'PL-5': 'GCATA', 'PM-5': 'GCATC', 'P4-5': 'GCATG', 'PN-5': 'GCATT', 'PO-5': 'GCCAA', 'P5-5': 'GCCAC',
               'PP-5': 'GCCAG', 'PQ-5': 'GCCAT', 'P6-5': 'GCCCA', 'PR-5': 'GCCCC', 'PS-5': 'GCCCG', 'P7-5': 'GCCCT',
               'PU-5': 'GCCGA', 'PV-5': 'GCCGC', 'P8-5': 'GCCGG', 'PW-5': 'GCCGT', 'PX-5': 'GCCTA', 'P9-5': 'GCCTC',
               'PY-5': 'GCCTG', 'PZ-5': 'GCCTT', 'QB-5': 'GCGAA', 'QD-5': 'GCGAC', 'Q0-5': 'GCGAG', 'QE-5': 'GCGAT',
               'QF-5': 'GCGCA', 'Q1-5': 'GCGCC', 'QH-5': 'GCGCG', 'QI-5': 'GCGCT', 'Q2-5': 'GCGGA', 'QJ-5': 'GCGGC',
               'QK-5': 'GCGGG', 'Q3-5': 'GCGGT', 'QL-5': 'GCGTA', 'QM-5': 'GCGTC', 'Q4-5': 'GCGTG', 'QN-5': 'GCGTT',
               'QO-5': 'GCTAA', 'Q5-5': 'GCTAC', 'QP-5': 'GCTAG', 'QQ-5': 'GCTAT', 'Q6-5': 'GCTCA', 'QR-5': 'GCTCC',
               'QS-5': 'GCTCG', 'Q7-5': 'GCTCT', 'QU-5': 'GCTGA', 'QV-5': 'GCTGC', 'Q8-5': 'GCTGG', 'QW-5': 'GCTGT',
               'QX-5': 'GCTTA', 'Q9-5': 'GCTTC', 'QY-5': 'GCTTG', 'QZ-5': 'GCTTT', '6B-5': 'GGAAA', '6D-5': 'GGAAC',
               '60-5': 'GGAAG', '6E-5': 'GGAAT', '6F-5': 'GGACA', '61-5': 'GGACC', '6H-5': 'GGACG', '6I-5': 'GGACT',
               '62-5': 'GGAGA', '6J-5': 'GGAGC', '6K-5': 'GGAGG', '63-5': 'GGAGT', '6L-5': 'GGATA', '6M-5': 'GGATC',
               '64-5': 'GGATG', '6N-5': 'GGATT', '6O-5': 'GGCAA', '65-5': 'GGCAC', '6P-5': 'GGCAG', '6Q-5': 'GGCAT',
               '66-5': 'GGCCA', '6R-5': 'GGCCC', '6S-5': 'GGCCG', '67-5': 'GGCCT', '6U-5': 'GGCGA', '6V-5': 'GGCGC',
               '68-5': 'GGCGG', '6W-5': 'GGCGT', '6X-5': 'GGCTA', '69-5': 'GGCTC', '6Y-5': 'GGCTG', '6Z-5': 'GGCTT',
               'RB-5': 'GGGAA', 'RD-5': 'GGGAC', 'R0-5': 'GGGAG', 'RE-5': 'GGGAT', 'RF-5': 'GGGCA', 'R1-5': 'GGGCC',
               'RH-5': 'GGGCG', 'RI-5': 'GGGCT', 'R2-5': 'GGGGA', 'RJ-5': 'GGGGC', 'RK-5': 'GGGGG', 'R3-5': 'GGGGT',
               'RL-5': 'GGGTA', 'RM-5': 'GGGTC', 'R4-5': 'GGGTG', 'RN-5': 'GGGTT', 'RO-5': 'GGTAA', 'R5-5': 'GGTAC',
               'RP-5': 'GGTAG', 'RQ-5': 'GGTAT', 'R6-5': 'GGTCA', 'RR-5': 'GGTCC', 'RS-5': 'GGTCG', 'R7-5': 'GGTCT',
               'RU-5': 'GGTGA', 'RV-5': 'GGTGC', 'R8-5': 'GGTGG', 'RW-5': 'GGTGT', 'RX-5': 'GGTTA', 'R9-5': 'GGTTC',
               'RY-5': 'GGTTG', 'RZ-5': 'GGTTT', 'SB-5': 'GTAAA', 'SD-5': 'GTAAC', 'S0-5': 'GTAAG', 'SE-5': 'GTAAT',
               'SF-5': 'GTACA', 'S1-5': 'GTACC', 'SH-5': 'GTACG', 'SI-5': 'GTACT', 'S2-5': 'GTAGA', 'SJ-5': 'GTAGC',
               'SK-5': 'GTAGG', 'S3-5': 'GTAGT', 'SL-5': 'GTATA', 'SM-5': 'GTATC', 'S4-5': 'GTATG', 'SN-5': 'GTATT',
               'SO-5': 'GTCAA', 'S5-5': 'GTCAC', 'SP-5': 'GTCAG', 'SQ-5': 'GTCAT', 'S6-5': 'GTCCA', 'SR-5': 'GTCCC',
               'SS-5': 'GTCCG', 'S7-5': 'GTCCT', 'SU-5': 'GTCGA', 'SV-5': 'GTCGC', 'S8-5': 'GTCGG', 'SW-5': 'GTCGT',
               'SX-5': 'GTCTA', 'S9-5': 'GTCTC', 'SY-5': 'GTCTG', 'SZ-5': 'GTCTT', '7B-5': 'GTGAA', '7D-5': 'GTGAC',
               '70-5': 'GTGAG', '7E-5': 'GTGAT', '7F-5': 'GTGCA', '71-5': 'GTGCC', '7H-5': 'GTGCG', '7I-5': 'GTGCT',
               '72-5': 'GTGGA', '7J-5': 'GTGGC', '7K-5': 'GTGGG', '73-5': 'GTGGT', '7L-5': 'GTGTA', '7M-5': 'GTGTC',
               '74-5': 'GTGTG', '7N-5': 'GTGTT', '7O-5': 'GTTAA', '75-5': 'GTTAC', '7P-5': 'GTTAG', '7Q-5': 'GTTAT',
               '76-5': 'GTTCA', '7R-5': 'GTTCC', '7S-5': 'GTTCG', '77-5': 'GTTCT', '7U-5': 'GTTGA', '7V-5': 'GTTGC',
               '78-5': 'GTTGG', '7W-5': 'GTTGT', '7X-5': 'GTTTA', '79-5': 'GTTTC', '7Y-5': 'GTTTG', '7Z-5': 'GTTTT',
               'UB-5': 'TAAAA', 'UD-5': 'TAAAC', 'U0-5': 'TAAAG', 'UE-5': 'TAAAT', 'UF-5': 'TAACA', 'U1-5': 'TAACC',
               'UH-5': 'TAACG', 'UI-5': 'TAACT', 'U2-5': 'TAAGA', 'UJ-5': 'TAAGC', 'UK-5': 'TAAGG', 'U3-5': 'TAAGT',
               'UL-5': 'TAATA', 'UM-5': 'TAATC', 'U4-5': 'TAATG', 'UN-5': 'TAATT', 'UO-5': 'TACAA', 'U5-5': 'TACAC',
               'UP-5': 'TACAG', 'UQ-5': 'TACAT', 'U6-5': 'TACCA', 'UR-5': 'TACCC', 'US-5': 'TACCG', 'U7-5': 'TACCT',
               'UU-5': 'TACGA', 'UV-5': 'TACGC', 'U8-5': 'TACGG', 'UW-5': 'TACGT', 'UX-5': 'TACTA', 'U9-5': 'TACTC',
               'UY-5': 'TACTG', 'UZ-5': 'TACTT', 'VB-5': 'TAGAA', 'VD-5': 'TAGAC', 'V0-5': 'TAGAG', 'VE-5': 'TAGAT',
               'VF-5': 'TAGCA', 'V1-5': 'TAGCC', 'VH-5': 'TAGCG', 'VI-5': 'TAGCT', 'V2-5': 'TAGGA', 'VJ-5': 'TAGGC',
               'VK-5': 'TAGGG', 'V3-5': 'TAGGT', 'VL-5': 'TAGTA', 'VM-5': 'TAGTC', 'V4-5': 'TAGTG', 'VN-5': 'TAGTT',
               'VO-5': 'TATAA', 'V5-5': 'TATAC', 'VP-5': 'TATAG', 'VQ-5': 'TATAT', 'V6-5': 'TATCA', 'VR-5': 'TATCC',
               'VS-5': 'TATCG', 'V7-5': 'TATCT', 'VU-5': 'TATGA', 'VV-5': 'TATGC', 'V8-5': 'TATGG', 'VW-5': 'TATGT',
               'VX-5': 'TATTA', 'V9-5': 'TATTC', 'VY-5': 'TATTG', 'VZ-5': 'TATTT', '8B-5': 'TCAAA', '8D-5': 'TCAAC',
               '80-5': 'TCAAG', '8E-5': 'TCAAT', '8F-5': 'TCACA', '81-5': 'TCACC', '8H-5': 'TCACG', '8I-5': 'TCACT',
               '82-5': 'TCAGA', '8J-5': 'TCAGC', '8K-5': 'TCAGG', '83-5': 'TCAGT', '8L-5': 'TCATA', '8M-5': 'TCATC',
               '84-5': 'TCATG', '8N-5': 'TCATT', '8O-5': 'TCCAA', '85-5': 'TCCAC', '8P-5': 'TCCAG', '8Q-5': 'TCCAT',
               '86-5': 'TCCCA', '8R-5': 'TCCCC', '8S-5': 'TCCCG', '87-5': 'TCCCT', '8U-5': 'TCCGA', '8V-5': 'TCCGC',
               '88-5': 'TCCGG', '8W-5': 'TCCGT', '8X-5': 'TCCTA', '89-5': 'TCCTC', '8Y-5': 'TCCTG', '8Z-5': 'TCCTT',
               'WB-5': 'TCGAA', 'WD-5': 'TCGAC', 'W0-5': 'TCGAG', 'WE-5': 'TCGAT', 'WF-5': 'TCGCA', 'W1-5': 'TCGCC',
               'WH-5': 'TCGCG', 'WI-5': 'TCGCT', 'W2-5': 'TCGGA', 'WJ-5': 'TCGGC', 'WK-5': 'TCGGG', 'W3-5': 'TCGGT',
               'WL-5': 'TCGTA', 'WM-5': 'TCGTC', 'W4-5': 'TCGTG', 'WN-5': 'TCGTT', 'WO-5': 'TCTAA', 'W5-5': 'TCTAC',
               'WP-5': 'TCTAG', 'WQ-5': 'TCTAT', 'W6-5': 'TCTCA', 'WR-5': 'TCTCC', 'WS-5': 'TCTCG', 'W7-5': 'TCTCT',
               'WU-5': 'TCTGA', 'WV-5': 'TCTGC', 'W8-5': 'TCTGG', 'WW-5': 'TCTGT', 'WX-5': 'TCTTA', 'W9-5': 'TCTTC',
               'WY-5': 'TCTTG', 'WZ-5': 'TCTTT', 'XB-5': 'TGAAA', 'XD-5': 'TGAAC', 'X0-5': 'TGAAG', 'XE-5': 'TGAAT',
               'XF-5': 'TGACA', 'X1-5': 'TGACC', 'XH-5': 'TGACG', 'XI-5': 'TGACT', 'X2-5': 'TGAGA', 'XJ-5': 'TGAGC',
               'XK-5': 'TGAGG', 'X3-5': 'TGAGT', 'XL-5': 'TGATA', 'XM-5': 'TGATC', 'X4-5': 'TGATG', 'XN-5': 'TGATT',
               'XO-5': 'TGCAA', 'X5-5': 'TGCAC', 'XP-5': 'TGCAG', 'XQ-5': 'TGCAT', 'X6-5': 'TGCCA', 'XR-5': 'TGCCC',
               'XS-5': 'TGCCG', 'X7-5': 'TGCCT', 'XU-5': 'TGCGA', 'XV-5': 'TGCGC', 'X8-5': 'TGCGG', 'XW-5': 'TGCGT',
               'XX-5': 'TGCTA', 'X9-5': 'TGCTC', 'XY-5': 'TGCTG', 'XZ-5': 'TGCTT', '9B-5': 'TGGAA', '9D-5': 'TGGAC',
               '90-5': 'TGGAG', '9E-5': 'TGGAT', '9F-5': 'TGGCA', '91-5': 'TGGCC', '9H-5': 'TGGCG', '9I-5': 'TGGCT',
               '92-5': 'TGGGA', '9J-5': 'TGGGC', '9K-5': 'TGGGG', '93-5': 'TGGGT', '9L-5': 'TGGTA', '9M-5': 'TGGTC',
               '94-5': 'TGGTG', '9N-5': 'TGGTT', '9O-5': 'TGTAA', '95-5': 'TGTAC', '9P-5': 'TGTAG', '9Q-5': 'TGTAT',
               '96-5': 'TGTCA', '9R-5': 'TGTCC', '9S-5': 'TGTCG', '97-5': 'TGTCT', '9U-5': 'TGTGA', '9V-5': 'TGTGC',
               '98-5': 'TGTGG', '9W-5': 'TGTGT', '9X-5': 'TGTTA', '99-5': 'TGTTC', '9Y-5': 'TGTTG', '9Z-5': 'TGTTT',
               'YB-5': 'TTAAA', 'YD-5': 'TTAAC', 'Y0-5': 'TTAAG', 'YE-5': 'TTAAT', 'YF-5': 'TTACA', 'Y1-5': 'TTACC',
               'YH-5': 'TTACG', 'YI-5': 'TTACT', 'Y2-5': 'TTAGA', 'YJ-5': 'TTAGC', 'YK-5': 'TTAGG', 'Y3-5': 'TTAGT',
               'YL-5': 'TTATA', 'YM-5': 'TTATC', 'Y4-5': 'TTATG', 'YN-5': 'TTATT', 'YO-5': 'TTCAA', 'Y5-5': 'TTCAC',
               'YP-5': 'TTCAG', 'YQ-5': 'TTCAT', 'Y6-5': 'TTCCA', 'YR-5': 'TTCCC', 'YS-5': 'TTCCG', 'Y7-5': 'TTCCT',
               'YU-5': 'TTCGA', 'YV-5': 'TTCGC', 'Y8-5': 'TTCGG', 'YW-5': 'TTCGT', 'YX-5': 'TTCTA', 'Y9-5': 'TTCTC',
               'YY-5': 'TTCTG', 'YZ-5': 'TTCTT', 'ZB-5': 'TTGAA', 'ZD-5': 'TTGAC', 'Z0-5': 'TTGAG', 'ZE-5': 'TTGAT',
               'ZF-5': 'TTGCA', 'Z1-5': 'TTGCC', 'ZH-5': 'TTGCG', 'ZI-5': 'TTGCT', 'Z2-5': 'TTGGA', 'ZJ-5': 'TTGGC',
               'ZK-5': 'TTGGG', 'Z3-5': 'TTGGT', 'ZL-5': 'TTGTA', 'ZM-5': 'TTGTC', 'Z4-5': 'TTGTG', 'ZN-5': 'TTGTT',
               'ZO-5': 'TTTAA', 'Z5-5': 'TTTAC', 'ZP-5': 'TTTAG', 'ZQ-5': 'TTTAT', 'Z6-5': 'TTTCA', 'ZR-5': 'TTTCC',
               'ZS-5': 'TTTCG', 'Z7-5': 'TTTCT', 'ZU-5': 'TTTGA', 'ZV-5': 'TTTGC', 'Z8-5': 'TTTGG', 'ZW-5': 'TTTGT',
               'ZX-5': 'TTTTA', 'Z9-5': 'TTTTC', 'ZY-5': 'TTTTG', 'ZZ-5': 'TTTTT', 'B-1': 'A', 'D-1': 'C',
               '0-1': 'G', 'E-1': 'T', 'F-2': 'AA', '1-2': 'AC', 'H-2': 'AG', 'I-2': 'AT',
               '2-2': 'CA', 'J-2': 'CC', 'K-2': 'CG', '3-2': 'CT', 'L-2': 'GA', 'M-2': 'GC',
               '4-2': 'GG', 'N-2': 'GT', 'O-2': 'TA', '5-2': 'TC', 'P-2': 'TG', 'Q-2': 'TT',
               '6-3': 'AAA', 'R-3': 'AAC', 'S-3': 'AAG', '7-3': 'AAT', 'U-3': 'ACA', 'V-3': 'ACC',
               '8-3': 'ACG', 'W-3': 'ACT', 'X-3': 'AGA', '9-3': 'AGC', 'Y-3': 'AGG', 'Z-3': 'AGT',
               'DB-3': 'ATA', 'DD-3': 'ATC', 'D0-3': 'ATG', 'DE-3': 'ATT', 'DF-3': 'CAA', 'D1-3': 'CAC',
               'DH-3': 'CAG', 'DI-3': 'CAT', 'D2-3': 'CCA', 'DJ-3': 'CCC', 'DK-3': 'CCG', 'D3-3': 'CCT',
               'DL-3': 'CGA', 'DM-3': 'CGC', 'D4-3': 'CGG', 'DN-3': 'CGT', 'DO-3': 'CTA', 'D5-3': 'CTC',
               'DP-3': 'CTG', 'DQ-3': 'CTT', 'D6-3': 'GAA', 'DR-3': 'GAC', 'DS-3': 'GAG', 'D7-3': 'GAT',
               'DU-3': 'GCA', 'DV-3': 'GCC', 'D8-3': 'GCG', 'DW-3': 'GCT', 'DX-3': 'GGA', 'D9-3': 'GGC',
               'DY-3': 'GGG', 'DZ-3': 'GGT', '0B-3': 'GTA', '0D-3': 'GTC', '00-3': 'GTG', '0E-3': 'GTT',
               '0F-3': 'TAA', '01-3': 'TAC', '0H-3': 'TAG', '0I-3': 'TAT', '02-3': 'TCA', '0J-3': 'TCC',
               '0K-3': 'TCG', '03-3': 'TCT', '0L-3': 'TGA', '0M-3': 'TGC', '04-3': 'TGG', '0N-3': 'TGT',
               '0O-3': 'TTA', '05-3': 'TTC', '0P-3': 'TTG', '0Q-3': 'TTT', '06-4': 'AAAA', '0R-4': 'AAAC',
               '0S-4': 'AAAG', '07-4': 'AAAT', '0U-4': 'AACA', '0V-4': 'AACC', '08-4': 'AACG', '0W-4': 'AACT',
               '0X-4': 'AAGA', '09-4': 'AAGC', '0Y-4': 'AAGG', '0Z-4': 'AAGT', 'EB-4': 'AATA', 'ED-4': 'AATC',
               'E0-4': 'AATG', 'EE-4': 'AATT', 'EF-4': 'ACAA', 'E1-4': 'ACAC', 'EH-4': 'ACAG', 'EI-4': 'ACAT',
               'E2-4': 'ACCA', 'EJ-4': 'ACCC', 'EK-4': 'ACCG', 'E3-4': 'ACCT', 'EL-4': 'ACGA', 'EM-4': 'ACGC',
               'E4-4': 'ACGG', 'EN-4': 'ACGT', 'EO-4': 'ACTA', 'E5-4': 'ACTC', 'EP-4': 'ACTG', 'EQ-4': 'ACTT',
               'E6-4': 'AGAA', 'ER-4': 'AGAC', 'ES-4': 'AGAG', 'E7-4': 'AGAT', 'EU-4': 'AGCA', 'EV-4': 'AGCC',
               'E8-4': 'AGCG', 'EW-4': 'AGCT', 'EX-4': 'AGGA', 'E9-4': 'AGGC', 'EY-4': 'AGGG', 'EZ-4': 'AGGT',
               'FB-4': 'AGTA', 'FD-4': 'AGTC', 'F0-4': 'AGTG', 'FE-4': 'AGTT', 'FF-4': 'ATAA', 'F1-4': 'ATAC',
               'FH-4': 'ATAG', 'FI-4': 'ATAT', 'F2-4': 'ATCA', 'FJ-4': 'ATCC', 'FK-4': 'ATCG', 'F3-4': 'ATCT',
               'FL-4': 'ATGA', 'FM-4': 'ATGC', 'F4-4': 'ATGG', 'FN-4': 'ATGT', 'FO-4': 'ATTA', 'F5-4': 'ATTC',
               'FP-4': 'ATTG', 'FQ-4': 'ATTT', 'F6-4': 'CAAA', 'FR-4': 'CAAC', 'FS-4': 'CAAG', 'F7-4': 'CAAT',
               'FU-4': 'CACA', 'FV-4': 'CACC', 'F8-4': 'CACG', 'FW-4': 'CACT', 'FX-4': 'CAGA', 'F9-4': 'CAGC',
               'FY-4': 'CAGG', 'FZ-4': 'CAGT', '1B-4': 'CATA', '1D-4': 'CATC', '10-4': 'CATG', '1E-4': 'CATT',
               '1F-4': 'CCAA', '11-4': 'CCAC', '1H-4': 'CCAG', '1I-4': 'CCAT', '12-4': 'CCCA', '1J-4': 'CCCC',
               '1K-4': 'CCCG', '13-4': 'CCCT', '1L-4': 'CCGA', '1M-4': 'CCGC', '14-4': 'CCGG', '1N-4': 'CCGT',
               '1O-4': 'CCTA', '15-4': 'CCTC', '1P-4': 'CCTG', '1Q-4': 'CCTT', '16-4': 'CGAA', '1R-4': 'CGAC',
               '1S-4': 'CGAG', '17-4': 'CGAT', '1U-4': 'CGCA', '1V-4': 'CGCC', '18-4': 'CGCG', '1W-4': 'CGCT',
               '1X-4': 'CGGA', '19-4': 'CGGC', '1Y-4': 'CGGG', '1Z-4': 'CGGT', 'HB-4': 'CGTA', 'HD-4': 'CGTC',
               'H0-4': 'CGTG', 'HE-4': 'CGTT', 'HF-4': 'CTAA', 'H1-4': 'CTAC', 'HH-4': 'CTAG', 'HI-4': 'CTAT',
               'H2-4': 'CTCA', 'HJ-4': 'CTCC', 'HK-4': 'CTCG', 'H3-4': 'CTCT', 'HL-4': 'CTGA', 'HM-4': 'CTGC',
               'H4-4': 'CTGG', 'HN-4': 'CTGT', 'HO-4': 'CTTA', 'H5-4': 'CTTC', 'HP-4': 'CTTG', 'HQ-4': 'CTTT',
               'H6-4': 'GAAA', 'HR-4': 'GAAC', 'HS-4': 'GAAG', 'H7-4': 'GAAT', 'HU-4': 'GACA', 'HV-4': 'GACC',
               'H8-4': 'GACG', 'HW-4': 'GACT', 'HX-4': 'GAGA', 'H9-4': 'GAGC', 'HY-4': 'GAGG', 'HZ-4': 'GAGT',
               'IB-4': 'GATA', 'ID-4': 'GATC', 'I0-4': 'GATG', 'IE-4': 'GATT', 'IF-4': 'GCAA', 'I1-4': 'GCAC',
               'IH-4': 'GCAG', 'II-4': 'GCAT', 'I2-4': 'GCCA', 'IJ-4': 'GCCC', 'IK-4': 'GCCG', 'I3-4': 'GCCT',
               'IL-4': 'GCGA', 'IM-4': 'GCGC', 'I4-4': 'GCGG', 'IN-4': 'GCGT', 'IO-4': 'GCTA', 'I5-4': 'GCTC',
               'IP-4': 'GCTG', 'IQ-4': 'GCTT', 'I6-4': 'GGAA', 'IR-4': 'GGAC', 'IS-4': 'GGAG', 'I7-4': 'GGAT',
               'IU-4': 'GGCA', 'IV-4': 'GGCC', 'I8-4': 'GGCG', 'IW-4': 'GGCT', 'IX-4': 'GGGA', 'I9-4': 'GGGC',
               'IY-4': 'GGGG', 'IZ-4': 'GGGT', '2B-4': 'GGTA', '2D-4': 'GGTC', '20-4': 'GGTG', '2E-4': 'GGTT',
               '2F-4': 'GTAA', '21-4': 'GTAC', '2H-4': 'GTAG', '2I-4': 'GTAT', '22-4': 'GTCA', '2J-4': 'GTCC',
               '2K-4': 'GTCG', '23-4': 'GTCT', '2L-4': 'GTGA', '2M-4': 'GTGC', '24-4': 'GTGG', '2N-4': 'GTGT',
               '2O-4': 'GTTA', '25-4': 'GTTC', '2P-4': 'GTTG', '2Q-4': 'GTTT', '26-4': 'TAAA', '2R-4': 'TAAC',
               '2S-4': 'TAAG', '27-4': 'TAAT', '2U-4': 'TACA', '2V-4': 'TACC', '28-4': 'TACG', '2W-4': 'TACT',
               '2X-4': 'TAGA', '29-4': 'TAGC', '2Y-4': 'TAGG', '2Z-4': 'TAGT', 'JB-4': 'TATA', 'JD-4': 'TATC',
               'J0-4': 'TATG', 'JE-4': 'TATT', 'JF-4': 'TCAA', 'J1-4': 'TCAC', 'JH-4': 'TCAG', 'JI-4': 'TCAT',
               'J2-4': 'TCCA', 'JJ-4': 'TCCC', 'JK-4': 'TCCG', 'J3-4': 'TCCT', 'JL-4': 'TCGA', 'JM-4': 'TCGC',
               'J4-4': 'TCGG', 'JN-4': 'TCGT', 'JO-4': 'TCTA', 'J5-4': 'TCTC', 'JP-4': 'TCTG', 'JQ-4': 'TCTT',
               'J6-4': 'TGAA', 'JR-4': 'TGAC', 'JS-4': 'TGAG', 'J7-4': 'TGAT', 'JU-4': 'TGCA', 'JV-4': 'TGCC',
               'J8-4': 'TGCG', 'JW-4': 'TGCT', 'JX-4': 'TGGA', 'J9-4': 'TGGC', 'JY-4': 'TGGG', 'JZ-4': 'TGGT',
               'KB-4': 'TGTA', 'KD-4': 'TGTC', 'K0-4': 'TGTG', 'KE-4': 'TGTT', 'KF-4': 'TTAA', 'K1-4': 'TTAC',
               'KH-4': 'TTAG', 'KI-4': 'TTAT', 'K2-4': 'TTCA', 'KJ-4': 'TTCC', 'KK-4': 'TTCG', 'K3-4': 'TTCT',
               'KL-4': 'TTGA', 'KM-4': 'TTGC', 'K4-4': 'TTGG', 'KN-4': 'TTGT', 'KO-4': 'TTTA', 'K5-4': 'TTTC',
               'KP-4': 'TTTG', 'KQ-4': 'TTTT'}


def is_sequence(sequence):
    """
    Checks if sequence is only made up of ATCG's
    :param sequence: Sequence to be checked
    :return: boolean
    """
    valid = re.compile('^[ATCG]+$', re.IGNORECASE)
    return valid.match(sequence)


def is_license_plate(length, code):
    """
    Checks if the license plate is in the proper format
    :param length: Encoded sequence length
    :param code: license plate for the sequence
    :return: Boolean value of if it is valid or not
    """
    if not length.isdigit():
        return False

    if any(c in ('A', 'T', 'G', 'C') or (not c.isalpha() and not c.isdigit()) for c in code):
        return False

    return True


def encode_sequence(sequence, prefix):
    """
    Encodes the sequence into its corresponding license plate with given prefix (if given one)
    :param sequence: The sequence being encoded
    :param prefix: The prefix to use for the license plate
    :return: The license plate it encodes to
    """
    length = len(sequence)
    # Encode label
    if prefix is '':
        final_result = [(str(length) + '-')]
    else:
        final_result = [prefix + "-" + str(length) + "-"]

    work_sequence = sequence
    while work_sequence != '':
        try:
            final_result.append(encode_hash[work_sequence[0:5]])
        except KeyError as err:
            if not err.args:
                err.args = ('',)
            err.args = err.args + ("Error, exiting: Segment '" + work_sequence[0:5] +
                                   "' from sequence '" + sequence + "' is invalid.",)
            raise
        work_sequence = work_sequence[5:]

    return ''.join(final_result)


def decode_sequence(plate):
    """
    Decode the license-plate using the lookup table
    :param plate: license plate being decoded
    :return: The sequence it decodes to
    """
    if plate.count('-') == 2:
        fields = plate.split('-', 2)
        # Discard prefix, it isn't used at all
        fields = fields[1:]

    elif plate.count('-') == 1:
        fields = plate.split('-', 1)
    else:
        raise KeyError("Error, exiting: Provided license plate '" + plate + "' is not in a valid format.")

    if not is_license_plate(fields[0], fields[1]):
        raise KeyError("Error, exiting: Provided license plate '" + plate + "' is not in a valid format.")

    length = int(fields[0])
    code = fields[1]

    # Retrieve sequence
    remainder = length
    raw_result = []
    while code != '' and remainder > 0:
        if remainder >= 5:
            try:
                raw_result.append(decode_hash[code[0:2] + '-5'])
            except KeyError as err:
                if not err.args:
                    err.args = ('',)
                err.args = err.args + ("Error, exiting: Segment '" + code[0:2] + "' from license plate '"
                                       + plate + "' is invalid.",)
                raise
        else:
            try:
                raw_result.append(decode_hash[code[0:2] + '-' + str(remainder)])
            except KeyError as err:
                if not err.args:
                    err.args = ('',)
                err.args = err.args + ("Error, exiting: Segment '" + code[0:2] + "' from license plate '"
                                       + plate + "' is invalid.",)
                raise

        remainder -= 5
        code = code[2:]

    # Check if label make sense
    final_result = ''.join(raw_result)
    if len(final_result) != length or code != '':
        raise KeyError("Error, exiting: Invalid license plate '" + plate + "'.")

    return final_result


def convert(seq, encode, prefix):
    if prefix is None:
        prefix = ''
    else:
        if '-' in prefix or ' ' in prefix:
            sys.stderr.write("Warning: Dashes and spaces are not permitted in the license plate prefix."
                             "Program will remove all instances automatically from prefix '" + prefix + "'.\n")
            prefix = prefix.replace('-', '')
            prefix = prefix.replace(' ', '')

    if encode:
        # Encode
        if seq != '':
            cleaned = seq.upper().replace('U', 'T')
            if is_sequence(cleaned):
                return encode_sequence(cleaned, prefix)
            else:
                raise KeyError('Error, exiting: Illegal characters in line "' + seq + '"')
    else:
        # Decode
        if seq != '':
            cleaned = seq.upper()
            return decode_sequence(cleaned)


def run_as_script():
    """
    Main function
    """
    from argparse import RawDescriptionHelpFormatter
    parser = argparse.ArgumentParser(
        description="Program to generate license plates from given sequence or decode license plate to a sequence. "
                    "An optional prefix can be provided to the program such as 'tRF', 'rRF', 'iso', etc, if needed.",
        formatter_class=RawDescriptionHelpFormatter)
    parser.add_argument('sequencefile', type=str, help='Sequences that are to be encoded or decoded')
    parser.add_argument('encode_choice', type=str, choices=['en', 'de'],
                        help='Choice of encoding (en) or decoding (de)')
    parser.add_argument('--p', '--prefix', type=str, default=None, dest='prefix',
                        help="Prefix to be used for license plates. 'tRF', 'rRF', 'iso', etc")

    args = parser.parse_args()

    encode = None

    if args.encode_choice == 'en':
        encode = True
    elif args.encode_choice == 'de':
        encode = False

    illegal_characters = ['\n', '\r', '\t', '\0', '\f', '`', '?', '*', '\\', '<', '>', '|', '\"', ': ']

    for item in illegal_characters:
        if item in args.sequencefile.split('/')[-1]:
            raise KeyError('Error, exiting: Illegal character "' + item + '" in filename')

    try:
        sequence_file = open(args.sequencefile, 'r')
    except IOError:
        sys.stderr.write('Error, exiting: File ' + args.sequencefile + ' not found.\n')
        sys.exit(1)

    if sequence_file is None:
        sys.stderr.write('Error: An unexpected error has occurred, please ensure that the file ' + args.sequencefile +
                         ' exists\n')
        sys.exit(1)

    sequences = []
    for line in sequence_file:
        sequences.append(line.rstrip())
    sequence_file.close()

    done = {}

    for seq in sequences:
        try:
            done[seq] = convert(seq, encode, args.prefix)
        except KeyError as error:
            print(error)
            sys.exit(1)

    for item in done:
        sys.stdout.write(item + '\t' + done[item] + '\n')


if __name__ == "__main__":
    run_as_script()
