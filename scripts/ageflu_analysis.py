#############################################
#  Python script for association analyses   #
#############################################

from __future__ import division

# parse fasta file
def parsefasta(file, check_alg=1, check_dna=1, check_codon=1):
    # check if file is of FASTA format
    fhandle = filter(None, open(file, 'rU').readlines())
    if not re.search('^>', fhandle[0]): sys.exit('\nERROR: Incorrect sequence file format.\n')

    result_dict = {}
    for key, group in itertools.groupby(fhandle, lambda _: re.search('^>',_)):
        if key:
            header = group.next().strip().replace('>','')
        else:
            sequence = ''.join(map(lambda _:_.strip(),list(group)))
            # dna sequences
            if check_dna == 1 and set(list(sequence))&set(list('rhkdesqpvilmfyw')):
                sys.exit('\nERROR: Input must be DNA FASTA.\n')
            # codon-alignment
            if check_codon == 1 and len(sequence)%3 > 0:
                sys.exit('\nERROR: DNA sequence file given must be codon-aligned.\n')
            result_dict[header] = sequence
    # alignment
    if check_alg == 1 and len(set(map(len, result_dict.values()))) != 1:
        sys.exit('\nERROR: Input sequence file must be an alignment.\n')
    return result_dict

# translate codon-aligned nucleotide sequences
def translatedDNA(dna):
    dnacodontable = {'ttt':'F', 'ttc':'F', 'tta':'L', 'ttg':'L', 'ctt':'L', 'ctc':'L', 'cta':'L', 'ctg':'L', 'att':'I', 'atc':'I', 'ata':'I', 'atg':'M', 'gtt':'V', 'gtc':'V', 'gta':'V', 'gtg':'V', 'tct':'S', 'tcc':'S', 'tca':'S', 'tcg':'S', 'cct':'P', 'ccc':'P', 'cca':'P', 'ccg':'P', 'act':'T', 'acc':'T', 'aca':'T', 'acg':'T', 'gct':'A', 'gcc':'A', 'gca':'A', 'gcg':'A', 'tat':'Y', 'tac':'Y', 'taa':'stop', 'tag':'stop', 'cat':'H', 'cac':'H', 'caa':'Q', 'cag':'Q', 'aat':'N', 'aac':'N', 'aaa':'K', 'aag':'K', 'gat':'D', 'gac':'D', 'gaa':'E', 'gag':'E', 'tgt':'C', 'tgc':'C', 'tga':'stop', 'tgg':'W', 'cgt':'R', 'cgc':'R', 'cga':'R', 'cgg':'R', 'agt':'S', 'agc':'S', 'aga':'R', 'agg':'R', 'ggt':'G', 'ggc':'G', 'gga':'G', 'ggg':'G'}

    protein = []
    for _ in xrange(0,len(dna),3):
        codon = dna[_:_+3]
        if codon in dnacodontable:
            if dnacodontable[codon] == 'stop': break
            protein.append(dnacodontable[codon])
        elif re.match('---',codon): protein.append('-')
        else: protein.append('X')
    return ''.join(protein)

# check for N-X-S/T motif changes
def check_glyco(positions, anc_seq, desc_seq):
    glyco_sites = []
    for pos in positions:
        if anc_seq[pos-1] == 'N' and re.search('N[^PX-][ST]', anc_seq[pos-1:pos+2]) and not re.search('N[^PX-][ST]', desc_seq[pos-1:pos+2]):
            glyco_sites.append(pos)
        elif re.search('[^PX-]', anc_seq[pos-1]) and re.search('N[^PX-][ST]', anc_seq[pos-2:pos+1]) and not re.search('N[^PX-][ST]', desc_seq[pos-2:pos+1]):
            glyco_sites.append(pos)
        elif re.search('[ST]', anc_seq[pos-1]) and re.search('N[^PX-][ST]', anc_seq[pos-3:pos]) and not re.search('N[^PX-][ST]', desc_seq[pos-3:pos]):
            glyco_sites.append(pos)

    return glyco_sites

# odds ratio analyses
def calculateOR(n11, n10, n01, n00):

    ORnum = n11*n00
    ORden = n01*n10

    if ORnum == ORden == 0:
        OR_uMLE, wald_MLE_CI95L, wald_MLE_CI95R = 'NaN', '', ''
        wald_MLE_CI90L, wald_MLE_CI90R = '', ''
    elif ORnum == 0:
        OR_uMLE, wald_MLE_CI95L, wald_MLE_CI95R = 0, '', ''
        wald_MLE_CI90L, wald_MLE_CI90R = '', ''
    elif ORden == 0:
        OR_uMLE, wald_MLE_CI95L, wald_MLE_CI95R = 'InF', '', ''
        wald_MLE_CI90L, wald_MLE_CI90R = '', ''
    else:
        OR_uMLE = ORnum/ORden
        # Wald test & CI
        lnOR = numpy.log(OR_uMLE)
        SE_lnOR = numpy.sqrt(sum(map(lambda _:1/_ ,[n11,n10,n01,n00])))
        wald_MLE_CI95L, wald_MLE_CI95R =  numpy.exp(lnOR - 1.96*SE_lnOR), numpy.exp(lnOR + 1.96*SE_lnOR)
        wald_MLE_CI90L, wald_MLE_CI90R =  numpy.exp(lnOR - 1.645*SE_lnOR), numpy.exp(lnOR + 1.645*SE_lnOR)

    if n11 + n10 + n01 + n00 >= 1000:
        # G-test
        cmd = ['Rscript', "ageflu_gtest.R"] + map(str,[n11, n10, n01, n00])
        Routput = subprocess.check_output(cmd, universal_newlines=True, stderr=subprocess.PIPE).split('\n')
        pval = float(Routput[1])
        # small sample size CI - NaN
        ss_CI90L, ss_CI90R, ss_CI95L, ss_CI95R = 'NaN', 'NaN', 'NaN', 'NaN'
    else:
        # Barnard's
        cmd = ['Rscript', "ageflu_barnard.R"] + map(str,[n11, n10, n01, n00])
        Routput = subprocess.check_output(cmd, universal_newlines=True, stderr=subprocess.PIPE).split('\n')
        pval = float(Routput[-2].replace('[1]','').split()[1])
        # Agresti & Min unconditional exact CI for small sample sizes
        cmd = ['Rscript', "agresti_and_min.R"] + map(str,[n11, n10, n01, n00, 0.9])
        Routput = filter(None, subprocess.check_output(cmd, universal_newlines=True, stderr=subprocess.PIPE).split('\n'))
        ss_CI90L, ss_CI90R = map(lambda _: _ if _ == 'Inf' else float(_), re.findall('(\d+\.\d+|\d+|Inf)', Routput[-1]))

        cmd = ['Rscript', "agresti_and_min.R"] + map(str,[n11, n10, n01, n00, 0.95])
        Routput = filter(None, subprocess.check_output(cmd, universal_newlines=True, stderr=subprocess.PIPE).split('\n'))
        ss_CI95L, ss_CI95R = map(lambda _: _ if _ == 'Inf' else float(_), re.findall('(\d+\.\d+|\d+|Inf)', Routput[-1]))

    return OR_uMLE, wald_MLE_CI90L, wald_MLE_CI90R, wald_MLE_CI95L, wald_MLE_CI95R, pval, ss_CI90L, ss_CI90R, ss_CI95L, ss_CI95R

# main
if __name__ == '__main__':
    import argparse
    params = argparse.ArgumentParser(description='Reads closest pairs output from ageflu_getpairs.py and perform age-stratified association analyses.')
    params.add_argument('-i', '--input_file', type=str, required=True, help='Closest pairs output file from ageflu_getpairs.py.')
    params.add_argument('-a', '--aln', type=str, required=True, help='Codon-aligned nucleotide fasta of sequences.')
    params.add_argument('-sp', '--same_passage', action='store_true', help='Analyse pairs with same passage histories only.')
    params = params.parse_args()

    # parameters
    # get query subtype
    import re
    try:
        query_subtype = re.search('(H3N2|pH1N1|H1N1pdm09|BVic|BYam)', params.input_file).group().upper()
    except:
        query_subtype = raw_input('\nWARNING: Can\'t parsed query subtype from tree name. Enter subtype (H3N2/pH1N1/BVic/BYam): ')
    if query_subtype == 'PH1N1':
        query_subtype = 'H1N1PDM09'

    # parameters to obtain reference positions
    import sys, itertools
    from os.path import expanduser
    ref_subtype_dictionary = {'H1N1PDM09':'H1', 'H3N2':'H3', 'BVIC':'B73', 'BYAM':'B73'}
    HA_ref_numbering_fdat = parsefasta('H1pdm09_H3_FluB_NumberingRef.fa', check_dna=0, check_codon=0)
    queryST_to_refST_to_AbNum_to_RefNum = {'H1N1PDM09': {'H1':{}, 'H3':{}, 'B73':{}}, 'H3N2': {'H1':{}, 'H3':{}, 'B73':{}}, 'BVIC': {'H1':{}, 'H3':{}, 'B73':{}}, 'BYAM': {'H1':{}, 'H3':{}, 'B73':{}}}
    queryST_to_RefStrain = {'H1N1PDM09':'A/Texas/04/2009(H1N1)_H1pdm09absolutenumbering', 'H3N2':'A/Aichi/2/1968(H3N2)_H3absolutenumbering', 'BVIC':'B/Brisbane/60/2008_BVicAbs', 'BYAM':'B/Phuket/3073/2013_BYamAbs'}
    RefStrain_to_RefSeq = {'H1':HA_ref_numbering_fdat['A/Texas/04/2009(H1N1)_H1pdm09absolutenumbering'], 'H3':HA_ref_numbering_fdat['A/Aichi/2/68(H3N2)_H3numbering'], 'B73':HA_ref_numbering_fdat['B/HK/8/73_BYamNumbering']}
    FluB_insertion_querypos_to_refpos = {177:162.1, 178:162.2, 179:162.3}

    for subtype in queryST_to_RefStrain.keys():
        for RefStrain, RefSeq in RefStrain_to_RefSeq.items():
            CurrRefPosNum, CurrQueryPosNum = 0, 0
            for _ in xrange(len(RefSeq)):
                RecordPosBIN = -1
                if RefSeq[_] != '-':
                    CurrRefPosNum += 1
                    RecordPosBIN += 1
                if HA_ref_numbering_fdat[queryST_to_RefStrain[subtype]][_] != '-':
                    CurrQueryPosNum += 1
                    RecordPosBIN += 1
                if RecordPosBIN > 0:
                    queryST_to_refST_to_AbNum_to_RefNum[subtype][RefStrain][CurrQueryPosNum] = CurrRefPosNum
                #162-163 insertion for Flu B (162.1, 162.2, 162.3)
                if (subtype == 'BVIC' or subtype == 'BYAM') and _ in range(177,180):
                    if HA_ref_numbering_fdat[queryST_to_RefStrain[subtype]][_] != '-':
                        queryST_to_refST_to_AbNum_to_RefNum[subtype][RefStrain][CurrQueryPosNum] = FluB_insertion_querypos_to_refpos[_]

    ### POSITION PARAMETERS - FluA (all based on H3 numbering)###
    # antigenic sites gathered from Stray, S. J., & Pittman, L. B. (2012). Subtype-and antigenic site-specific differences in biophysical influences on evolution of influenza virus hemagglutinin. Virology journal, 9(1), 1.
    H3ant = {"A": range(121,130) + range(131,139) + range(142,147)  + [140],
             "B": range(155,161) + range(188,191) + range(192,195) + range(196,200) + [186, 246, 247],
             "C": [49, 50, 53, 54, 271, 273, 275, 276, 278],
             "D": range(201,208) + range(216, 221) + range(225, 228) + [167, 214, 222, 223, 242],
             "E": range(79, 84) + [62, 63, 75, 78, 91, 92]}

    H1ant = {"Sa": range(165,168) + [128, 129, 163, 247, 248],
             "Sb": [156, 159, 189, 190, 192, 193, 196, 197, 198],
             "Ca1": range(240, 246) + [169, 173, 207, 212],
             "Ca2": [132, 140, 143, 144, 145, 149, 224, 225],
             "Cb": range(259, 266) + [56, 79, 80, 82, 83, 117, 149, 255, 256],
             "H1C": range(271, 277) + range(278, 282) + [90, 285]}

    # H3 rbs gathered from Yang, H., Carney, P. J., Chang, J. C., Guo, Z., Villanueva, J. M., & Stevens, J. (2015). Structure and receptor binding preferences of recombinant human A (H3N2) virus hemagglutinins. Virology, 477, 18-31 - based on H3 numbering
    H3RBS = [98, 131, 132, 133, 134, 135, 136, 137, 138, 139, 140, 144, 145, 146, 153, 154, 155, 156, 157, 158, 159, 183, 184, 185, 186, 187, 188, 189, 190, 191, 192, 193, 194, 195, 196] + range(219,229)

    #passage sites - based on H3 numbering
    weiwei22pas = [111, 126, 137, 138, 144, 145, 155, 156, 158, 159, 185, 186, 193, 194, 199, 219, 226, 229, 246, 248, 276, 290]
    RaphaelPas = [13, 96, 137, 156, 183, 186, 189, 190, 193, 194, 196, 203, 219, 222, 246, 248, 328, 332, 335, 337, 352, 353, 376, 447, 465, 504, 515]
    SebastianPas = [101, 103, 123, 129, 131, 133, 137, 140, 142, 144, 145, 156, 157, 158, 159, 165, 187, 188, 192, 193, 196, 197, 198, 218, 222, 225, 226, 227, 252, 262, 269, 273, 279, 312, 474, 476, 498, 45, 510, 529, 48]

    ### POSITION PARAMETERS - Flu B (all based on B/HK/8/73 numbering)###
    FluB_insertion_querypos_to_refpos = {177:162.1, 178:162.2, 179:162.3}
    FluBant = {"120loop":range(51,55) + range(56,60)+range(73,76)+range(121,124)+range(179,182)+range(283,286)+[129,137],
               "150loop":range(139,143) + range(145,152),
               "160loop":range(162,166) + range(167,172),
               "190helix":range(194,201) + range(235,241)+[205],
               "162-163insertion":[162.1, 162.2, 162.3]}
    FluBRBS = range(193,203) + range(237,243) + range(136,144) + [95, 158, 191, 202]


    if query_subtype == 'H3N2':
        AS_dict, RBS_list = H3ant, H3RBS
    elif query_subtype == 'H1N1PDM09':
        AS_dict, RBS_list = H1ant, H3RBS
    else:
        AS_dict, RBS_list = FluBant, FluBRBS
    # AS + RBS
    functional_sites = list(set(RBS_list)|set([k for v in AS_dict.values() for k in v]))

    # read fasta file
    fdat_nuc = parsefasta(params.aln)
    fdat_aa = {k:translatedDNA(v) for k,v in fdat_nuc.items()}
    header_to_passage = {header:re.search('(ORI|MDCK|SIAT|EGG)', header).group() for header in fdat_nuc.keys()}

    # read input file and count pairs/substitutions
    print ('...counting pairs and sorting substitutions...')
    pt_count, pt_to_AAsublen, pt_to_RBASsub, abspos_to_pt_to_count, abspos_to_mutation_to_pt_to_count, abspos_to_refpos = {}, {}, {}, {}, {}, {}
    all_potential_glycosylation_sites = []

    fhandle = filter(None, open(params.input_file, 'rU'))
    fhandle.pop(0)
    for line in fhandle:
        try:
            pt, anc, desc, substitutions, substitutions_ref_positions = line.strip().split('\t')
            substitutions = substitutions.split(',')
        except:
            pt, anc, desc = line.strip().split('\t')
            substitutions = []

        # same passage filter
        if params.same_passage and header_to_passage[anc] != header_to_passage[desc]:
            continue

        # add to pair type count
        try:
            pt_count[pt] += 1
        except:
            pt_count[pt] = 1

        # add to pair type substitution count
        try:
            pt_to_AAsublen[pt][len(substitutions)] += 1
        except:
            try:
                pt_to_AAsublen[pt][len(substitutions)] = 1
            except:
                pt_to_AAsublen[pt] = {len(substitutions):1}

        # continue if NC pairs
        if pt == 'NC':
            continue

        if len(substitutions) > 0:
            # check if any substitutions are potential glycosylation sites
            substitutions_abs_positions = [int(re.search('\d+', _).group()) for _ in substitutions]
            potential_glycosylation_sites = check_glyco(substitutions_abs_positions, fdat_aa[anc], fdat_aa[desc])
            all_potential_glycosylation_sites = list(set(all_potential_glycosylation_sites)|set(potential_glycosylation_sites))


            for p, abs_pos in enumerate(substitutions_abs_positions):
                # add to position substitution count
                try:
                    abspos_to_pt_to_count[abs_pos][pt] += 1
                except:
                    try:
                        abspos_to_pt_to_count[abs_pos][pt] = 1
                    except:
                        abspos_to_pt_to_count[abs_pos] = {pt:1}

                # add to pos-mutation dictionary
                substitution = substitutions[p]
                try:
                    abspos_to_mutation_to_pt_to_count[abs_pos][substitution][pt] += 1
                except:
                    try:
                        abspos_to_mutation_to_pt_to_count[abs_pos][substitution][pt] = 1
                    except:
                        try:
                            abspos_to_mutation_to_pt_to_count[abs_pos][substitution] = {pt:1}
                        except:
                            abspos_to_mutation_to_pt_to_count[abs_pos] = {substitution:{pt:1}}

            # count if any sites are either RBS/AS
            ref_positions = []
            for _ in substitutions_ref_positions.split(','):
                if _ == '-':
                    ref_positions.append('-')
                elif re.search('162(\.1|\.2|\.3)', _):
                    ref_positions.append(float(_))
                else:
                    ref_positions.append(int(_))
            for _, refpos in enumerate(ref_positions):
                abspos_to_refpos[substitutions_abs_positions[_]] = refpos

            if set(functional_sites)&set(ref_positions):
                try:
                    pt_to_RBASsub[pt] += 1
                except:
                    pt_to_RBASsub[pt] = 1

    from scipy import stats
    from math import sqrt
    import statsmodels.api as sm
    import subprocess, numpy
    # distribution analyses
    print ('...pair distribution analyses...')
    maxmut = int(re.search('_MM(\d+)_', params.input_file).group(1))
    outfname = '{}_SP{}.txt'.format(re.sub('/*[^/]+/', '', params.input_file.replace('ageflu_evol-closest-pairs.', '').replace('.txt', '')), '1' if params.same_passage else '0')

    with open('ageflu_pairs-distribution.{}'.format(outfname), 'w') as output:
        # pair distribution - frequencies
        output.write('Frequency distribution of pair types stratified by number of amino acid substitutions\nPT/Sub_#\tCA\tAC\tCC\tAA\tNC\tSub#_TOT\n')
        col_counts = [0]*(maxmut)
        sublen_to_row_counts = {}
        for sub_len in xrange(maxmut+1):
            try:
                CA_count = pt_to_AAsublen['CA'][sub_len]
            except:
                CA_count = 0
                pt_to_AAsublen['CA'][sub_len] = 0
            try:
                AC_count = pt_to_AAsublen['AC'][sub_len]
            except:
                AC_count = 0
                pt_to_AAsublen['AC'][sub_len] = 0
            try:
                CC_count = pt_to_AAsublen['CC'][sub_len]
            except:
                CC_count = 0
                pt_to_AAsublen['CC'][sub_len] = 0
            try:
                AA_count = pt_to_AAsublen['AA'][sub_len]
            except:
                AA_count = 0
                pt_to_AAsublen['AA'][sub_len] = 0
            try:
                NC_count = pt_to_AAsublen['NC'][sub_len]
            except:
                NC_count = 0
                pt_to_AAsublen['NC'][sub_len] = 0

            row_counts = [CA_count, AC_count, CC_count, AA_count, NC_count]
            sublen_to_row_counts[sub_len] = row_counts

            for k, count in enumerate(row_counts):
                col_counts[k] += count
            output.write('{}\t{}\t{}\n'.format(sub_len, '\t'.join(map(str, row_counts)), sum(row_counts)))
        output.write('PT_TOT\t{}\t\n\n'.format('\t'.join(map(str, col_counts))))

        # pairs distribution - proportions
        output.write('Distribution of pair types stratified by number of amino acid substitutions\nPT/Sub_#\tCA\tAC\tCC\tAA\tNC\n')
        for sub_len in xrange(maxmut+1):
            output.write('{}\t{}\n'.format(sub_len, '\t'.join(map(lambda _: '{}'.format(_), [sublen_to_row_counts[sub_len][k]/col_tot for k, col_tot in enumerate(col_counts)]))))
        output.write('\n')

        # Wilcoxon signed-rank tests between pair types
        output.write('Wilcoxon signed-rank test between pair types\nPT_I\tPT_J\tp-value\n')
        for pt_i, pt_j in itertools.combinations(pt_to_AAsublen.keys(), 2):
            try:
                pt_i_counts = [1000*(pt_to_AAsublen[pt_i][sub_len]/pt_count[pt_i]) for sub_len in xrange(maxmut+1)]
            except:
                pt_count[pt_i] = 0
                pt_i_counts = [0]*(maxmut+1)

            try:
                pt_j_counts = [1000*(pt_to_AAsublen[pt_j][sub_len]/pt_count[pt_j]) for sub_len in xrange(maxmut+1)]
            except:
                pt_count[pt_j] = 0
                pt_j_counts = [0]*(maxmut+1)
            output.write('{}\t{}\t{}\n'.format(pt_i, pt_j, stats.wilcoxon(pt_i_counts, pt_j_counts).pvalue))

        # end-in-adult v. end-in-child
        pt_i_counts = [1000*(pt_to_AAsublen['AA'][sub_len] + pt_to_AAsublen['CA'][sub_len])/(pt_count['AA'] + pt_count['CA']) for sub_len in xrange(maxmut+1)]
        pt_j_counts = [1000*(pt_to_AAsublen['CC'][sub_len] + pt_to_AAsublen['AC'][sub_len])/(pt_count['CC'] + pt_count['AC']) for sub_len in xrange(maxmut+1)]
        output.write('end-in-adult\tend-in-child\t{}\n\n'.format(stats.wilcoxon(pt_i_counts, pt_j_counts).pvalue))

        # Association analyses between pair types (end-in-adult v. end-in-child) and propensity for substitution(s) - whole HA
        n1t = pt_count['CA'] + pt_count['AA']
        n0t = pt_count['CC'] + pt_count['AC']

        n11 = sum([pt_to_AAsublen['CA'][_] for _ in pt_to_AAsublen['CA'] if _ > 0]) + sum([pt_to_AAsublen['AA'][_] for _ in pt_to_AAsublen['AA'] if _ > 0])
        n01 = sum([pt_to_AAsublen['CC'][_] for _ in pt_to_AAsublen['CC'] if _ > 0]) + sum([pt_to_AAsublen['AC'][_] for _ in pt_to_AAsublen['AC'] if _ > 0])

        output.write('Association analyses between pair-type (end-in-adult v. end-in-child) and propensity for substitution(s) - whole HA\nPT/SF\tSubstitution\tNo Substitution\nEnd-in-adult\t{}\t{}\nEnd-in-child\t{}\t{}\n'.format(n11, n1t-n11, n01, n0t-n01))
        OR_uMLE, wald_MLE_CI90L, wald_MLE_CI90R, wald_MLE_CI95L, wald_MLE_CI95R, pval, ss_CI90L, ss_CI90R, ss_CI95L, ss_CI95R = calculateOR(n11, n1t-n11, n01, n0t-n01)
        output.write('OR\t{}\tp-value\t{}\nWald_95CI\t{}\nWald_90CI\t{}\nA&M_90CI\t{}\nA&M_95CI\t{}\n\n'.format(OR_uMLE, pval, ','.join(map(lambda _: '{}'.format(_), [wald_MLE_CI90L, wald_MLE_CI90R])), ','.join(map(lambda _: '{}'.format(_), [wald_MLE_CI95L, wald_MLE_CI95R])), ','.join(map(lambda _: '{}'.format(_), [ss_CI90L, ss_CI90R])), ','.join(map(lambda _: '{}'.format(_), [ss_CI95L, ss_CI95R]))))

        # Association analyses between pair types (end-in-adult v. end-in-child) and propensity for substitution(s) - RBS/AS only
        # n1t and n0t stays the same
        try:
            n11 = pt_to_RBASsub['CA']
        except:
            n11 = 0

        try:
            n11 += pt_to_RBASsub['AA']
        except:
            pass

        try:
            n01 = pt_to_RBASsub['CC']
        except:
            n01 = 0

        try:
            n01 += pt_to_RBASsub['AC']
        except:
            pass

        output.write('Association analyses between pair-type (end-in-adult v. end-in-child) and propensity for substitution(s) - RBS/AS\nPT/SF\tSubstitution\tNo Substitution\nEnd-in-adult\t{}\t{}\nEnd-in-child\t{}\t{}\n'.format(n11, n1t-n11, n01, n0t-n01))
        OR_uMLE, wald_MLE_CI90L, wald_MLE_CI90R, wald_MLE_CI95L, wald_MLE_CI95R, pval, ss_CI90L, ss_CI90R, ss_CI95L, ss_CI95R = calculateOR(n11, n1t-n11, n01, n0t-n01)
        output.write('OR\t{}\tp-value\t{}\nWald_95CI\t{}\nWald_90CI\t{}\nA&M_90CI\t{}\nA&M_95CI\t{}\n\n'.format(OR_uMLE, pval, ','.join(map(lambda _: '{}'.format(_), [wald_MLE_CI90L, wald_MLE_CI90R])), ','.join(map(lambda _: '{}'.format(_), [wald_MLE_CI95L, wald_MLE_CI95R])), ','.join(map(lambda _: '{}'.format(_), [ss_CI90L, ss_CI90R])), ','.join(map(lambda _: '{}'.format(_), [ss_CI95L, ss_CI95R]))))

        # Association analyses between pair types (CA v. AC) and propensity for substitution(s) - whole HA
        n1t = pt_count['CA']
        n0t = pt_count['AC']

        n11 = sum([pt_to_AAsublen['CA'][_] for _ in pt_to_AAsublen['CA'] if _ > 0])
        n01 = sum([pt_to_AAsublen['AC'][_] for _ in pt_to_AAsublen['AC'] if _ > 0])

        output.write('Association analyses between pair-type (CA v. AC) and propensity for substitution(s) - whole HA\nPT/SF\tSubstitution\tNo Substitution\nCA\t{}\t{}\nAC\t{}\t{}\n'.format(n11, n1t-n11, n01, n0t-n01))
        OR_uMLE, wald_MLE_CI90L, wald_MLE_CI90R, wald_MLE_CI95L, wald_MLE_CI95R, pval, ss_CI90L, ss_CI90R, ss_CI95L, ss_CI95R = calculateOR(n11, n1t-n11, n01, n0t-n01)
        output.write('OR\t{}\tp-value\t{}\nWald_95CI\t{}\nWald_90CI\t{}\nA&M_90CI\t{}\nA&M_95CI\t{}\n\n'.format(OR_uMLE, pval, ','.join(map(lambda _: '{}'.format(_), [wald_MLE_CI90L, wald_MLE_CI90R])), ','.join(map(lambda _: '{}'.format(_), [wald_MLE_CI95L, wald_MLE_CI95R])), ','.join(map(lambda _: '{}'.format(_), [ss_CI90L, ss_CI90R])), ','.join(map(lambda _: '{}'.format(_), [ss_CI95L, ss_CI95R]))))

        # Association analyses between pair types (CA v. AC) and propensity for substitution(s) - RBS/AS only
        # n1t and n0t stays the same
        try:
            n11 = pt_to_RBASsub['CA']
        except:
            n11 = 0

        try:
            n01 = pt_to_RBASsub['AC']
        except:
            n01 = 0

        output.write('Association analyses between pair-type (CA v. AC) and propensity for substitution(s) - RBS/AS\nPT/SF\tSubstitution\tNo Substitution\nCA\t{}\t{}\nAC\t{}\t{}\n'.format(n11, n1t-n11, n01, n0t-n01))
        OR_uMLE, wald_MLE_CI90L, wald_MLE_CI90R, wald_MLE_CI95L, wald_MLE_CI95R, pval, ss_CI90L, ss_CI90R, ss_CI95L, ss_CI95R = calculateOR(n11, n1t-n11, n01, n0t-n01)
        output.write('OR\t{}\tp-value\t{}\nWald_95CI\t{}\nWald_90CI\t{}\nA&M_90CI\t{}\nA&M_95CI\t{}\n\n'.format(OR_uMLE, pval, ','.join(map(lambda _: '{}'.format(_), [wald_MLE_CI90L, wald_MLE_CI90R])), ','.join(map(lambda _: '{}'.format(_), [wald_MLE_CI95L, wald_MLE_CI95R])), ','.join(map(lambda _: '{}'.format(_), [ss_CI90L, ss_CI90R])), ','.join(map(lambda _: '{}'.format(_), [ss_CI95L, ss_CI95R]))))

    # Age-associated sites
    print ('...calculating age-associated sites...')
    abspos_to_OR_results = {}
    for abs_pos in abspos_to_pt_to_count.keys():
        # end-in-adult v. end-in-child
        try:
            n11 = abspos_to_pt_to_count[abs_pos]['CA']
            CA_11 = abspos_to_pt_to_count[abs_pos]['CA']
        except:
            n11 = 0
            CA_11 = 0

        try:
            n11 += abspos_to_pt_to_count[abs_pos]['AA']
        except:
            pass

        try:
            n01 = abspos_to_pt_to_count[abs_pos]['CC']
        except:
            n01 = 0

        try:
            n01 += abspos_to_pt_to_count[abs_pos]['AC']
            AC_11 = abspos_to_pt_to_count[abs_pos]['AC']
        except:
            AC_11 = 0

        n10 = pt_count['CA'] + pt_count['AA'] - n11
        n00 = pt_count['CC'] + pt_count['AC'] - n01

        abspos_to_OR_results[abs_pos] = {'all': calculateOR(n11, n10, n01, n00), 'caac': calculateOR(CA_11, pt_count['CA']-CA_11, AC_11, pt_count['AC']-AC_11)}

    # get q values
    all_OR_pvals = [abspos_to_OR_results[abs_pos]['all'][5] for abs_pos in abspos_to_OR_results.keys()]
    all_OR_qvals = sm.stats.multipletests(all_OR_pvals, method='fdr_bh')[1].tolist()
    all_pval_to_qval = {all_OR_pvals[_]:qvalue for _, qvalue in enumerate(all_OR_qvals)}

    caac_OR_pvals = [abspos_to_OR_results[abs_pos]['caac'][5] for abs_pos in abspos_to_OR_results.keys()]
    caac_OR_qvals = sm.stats.multipletests(caac_OR_pvals, method='fdr_bh')[1].tolist()
    caac_pval_to_qval = {caac_OR_pvals[_]:qvalue for _, qvalue in enumerate(caac_OR_qvals)}

    with open('ageflu_age-associated-sites.{}'.format(outfname), 'w') as output:
        output.write('abs_pos\tref(H3)\tref(H1)\tref(B73)\tRBS\tAS\tGS\tPS\t'
                     'OR_uMLE_(EIA v. EIC)\twald_MLE_CI90L\twald_MLE_CI90R\twald_MLE_CI95L\twald_MLE_CI95R\t'
                     'p-value_(EIA v. EIC)\tss_CI90L\tss_CI90R\tss_CI95L\tss_CI95R\tq-value_(EIA v. EIC)\t'
                     'OR_uMLE_(CA v. AC)\twald_MLE_CI90L\twald_MLE_CI90R\twald_MLE_CI95L\twald_MLE_CI95R\t'
                     'p-value_(CA v. AC)\tss_CI90L\tss_CI90R\tss_CI95L\tss_CI95R\tq-value_(CA v. AC)\t'
                     'Mutations({})\tCA\tAC\tAA\tCC\tend-in-adult\tend-in-child\tTotal\n'.format(ref_subtype_dictionary[query_subtype]))
        for abs_pos in abspos_to_OR_results.keys():
            # reference positions
            try:
                main_ref_pos = queryST_to_refST_to_AbNum_to_RefNum[query_subtype][ref_subtype_dictionary[query_subtype]][abs_pos]
            except:
                main_ref_pos = '-'
            try:
                H3_ref_pos = queryST_to_refST_to_AbNum_to_RefNum[query_subtype]['H3'][abs_pos]
            except:
                H3_ref_pos = '-'
            try:
                H1_ref_pos = queryST_to_refST_to_AbNum_to_RefNum[query_subtype]['H1'][abs_pos]
            except:
                H1_ref_pos = '-'
            try:
                B_ref_pos = queryST_to_refST_to_AbNum_to_RefNum[query_subtype]['B73'][abs_pos]
            except:
                B_ref_pos = '-'

            if re.search('(H1N1PDM09|H3N2)', query_subtype):
                RBS_binary = 1 if H3_ref_pos in RBS_list else ''
            else:
                RBS_binary = 1 if main_ref_pos in RBS_list else ''

            try:
                if re.search('(H1N1PDM09|H3N2)', query_subtype):
                    AS_type = [_ for _ in AS_dict if H3_ref_pos in AS_dict[_]][0]
                else:
                    AS_type = [_ for _ in AS_dict if main_ref_pos in AS_dict[_]][0]
            except:
                AS_type = ''
            GS_binary = 1 if abs_pos in all_potential_glycosylation_sites else ''
            PS_binary = 1 if H3_ref_pos in SebastianPas else ''

            # OR results
            all_OR_result = map(str, ['{:4f}'.format(_) if isinstance(_, float) else _ for _ in abspos_to_OR_results[abs_pos]['all']])
            caac_OR_result = map(str, ['{:4f}'.format(_) if isinstance(_, float) else _ for _ in abspos_to_OR_results[abs_pos]['caac']])

            # substitutions
            substitution_line = ['{}({})'.format(re.sub(str(abs_pos), str(main_ref_pos), substitution), ','.join(['{}:{}'.format(pt, count) for pt, count in abspos_to_mutation_to_pt_to_count[abs_pos][substitution].items()])) for substitution in abspos_to_mutation_to_pt_to_count[abs_pos].keys()]

            # pair counts
            count_line = []
            end_in_adult_count, end_in_child_count = 0, 0
            try:
                count_line.append(abspos_to_pt_to_count[abs_pos]['CA'])
                end_in_adult_count += abspos_to_pt_to_count[abs_pos]['CA']
            except:
                count_line.append(0)

            try:
                count_line.append(abspos_to_pt_to_count[abs_pos]['AC'])
                end_in_child_count += abspos_to_pt_to_count[abs_pos]['AC']
            except:
                count_line.append(0)

            try:
                count_line.append(abspos_to_pt_to_count[abs_pos]['AA'])
                end_in_adult_count += abspos_to_pt_to_count[abs_pos]['AA']
            except:
                count_line.append(0)

            try:
                count_line.append(abspos_to_pt_to_count[abs_pos]['CC'])
                end_in_child_count += abspos_to_pt_to_count[abs_pos]['CC']
            except:
                count_line.append(0)

            count_line = count_line + [end_in_adult_count, end_in_child_count, end_in_adult_count+end_in_child_count]

            write_line = map(str, [abs_pos, H3_ref_pos, H1_ref_pos, B_ref_pos, RBS_binary, AS_type, GS_binary, PS_binary, '\t'.join(all_OR_result), all_pval_to_qval[abspos_to_OR_results[abs_pos]['all'][5]], '\t'.join(caac_OR_result), caac_pval_to_qval[abspos_to_OR_results[abs_pos]['caac'][5]], ';'.join(substitution_line), '\t'.join(map(str, count_line))])

            output.write('{}\n'.format('\t'.join(write_line)))
            output.flush()

    print ('...done.\n')