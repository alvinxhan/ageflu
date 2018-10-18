from __future__ import division
from pathos.helpers import mp
from os.path import expanduser as eu
from scipy import stats
import re
import sys
import argparse
import itertools
import subprocess
import statsmodels.api as sm
import time
import numpy

def get_contingency_counts(inputdict):
    # end-in-adult x 1
    try:
        p11 = inputdict['CA']['1']
    except:
        p11 = 0

    try:
        p11 += inputdict['AA']['1']
    except:
        pass

    # end-in-child x 1
    try:
        p01 = inputdict['CC']['1']
    except:
        p01 = 0

    try:
        p01 += inputdict['AC']['1']
    except:
        pass

    # end-in-adult x 0
    try:
        p10 = inputdict['CA']['0']
    except:
        p10 = 0

    try:
        p10 += inputdict['AA']['0']
    except:
        pass

    # end-in-child x 0
    try:
        p00 = inputdict['CC']['0']
    except:
        p00 = 0

    try:
        p00 += inputdict['AC']['0']
    except:
        pass

    return p11, p10, p01, p00

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
        # Wald test\tCI
        lnOR = numpy.log(OR_uMLE)
        SE_lnOR = numpy.sqrt(sum(map(lambda _:1/_ ,[n11,n10,n01,n00])))
        wald_MLE_CI95L, wald_MLE_CI95R =  numpy.exp(lnOR - 1.96*SE_lnOR), numpy.exp(lnOR + 1.96*SE_lnOR)
        wald_MLE_CI90L, wald_MLE_CI90R =  numpy.exp(lnOR - 1.645*SE_lnOR), numpy.exp(lnOR + 1.645*SE_lnOR)

    if n11 + n10 + n01 + n00 >= 1000:
        # G-test
        cmd = ['Rscript', eu("~/Dropbox/age_final_submission/ageflu/scripts/ageflu_gtest.R")] + map(str,[n11, n10, n01, n00])
        Routput = subprocess.check_output(cmd, universal_newlines=True, stderr=subprocess.PIPE).split('\n')
        pval = float(Routput[1])
        # small sample size CI - NaN
        ss_CI90L, ss_CI90R, ss_CI95L, ss_CI95R = 'NaN', 'NaN', 'NaN', 'NaN'
    else:
        # Barnard's
        cmd = ['Rscript', eu("~/Dropbox/age_final_submission/ageflu/scripts/ageflu_barnard.R")] + map(str,[n11, n10, n01, n00])
        Routput = subprocess.check_output(cmd, universal_newlines=True, stderr=subprocess.PIPE).split('\n')
        pval = float(Routput[-2].replace('[1]','').split()[1])
        # Agresti\tMin unconditional exact CI for small sample sizes
        cmd = ['Rscript', eu("~/Dropbox/age_final_submission/ageflu/scripts/agresti_and_min.R")] + map(str,[n11, n10, n01, n00, 0.9])
        Routput = filter(None, subprocess.check_output(cmd, universal_newlines=True, stderr=subprocess.PIPE).split('\n'))
        ss_CI90L, ss_CI90R = map(lambda _: _ if _ == 'Inf' else float(_), re.findall('(\d+\.\d+|\d+|Inf)', Routput[-1]))

        cmd = ['Rscript', eu("~/Dropbox/age_final_submission/ageflu/scripts/agresti_and_min.R")] + map(str,[n11, n10, n01, n00, 0.95])
        Routput = filter(None, subprocess.check_output(cmd, universal_newlines=True, stderr=subprocess.PIPE).split('\n'))
        ss_CI95L, ss_CI95R = map(lambda _: _ if _ == 'Inf' else float(_), re.findall('(\d+\.\d+|\d+|Inf)', Routput[-1]))

    return OR_uMLE, wald_MLE_CI90L, wald_MLE_CI90R, wald_MLE_CI95L, wald_MLE_CI95R, pval, ss_CI90L, ss_CI90R, ss_CI95L, ss_CI95R


###---FUNCTION: PARSE DNA CODON-ALIGNED FASTA---###
def parsefasta(file, check_alg=1, check_dna=1, check_codon=1):
    # check if file is of FASTA format
    fhandle = filter(None, open(file, 'rU').readlines())
    if not re.search('^>', fhandle[0]): sys.exit('\nERROR: Incorrect sequence file format.\n')

    result = {}
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
            result[header] = sequence
    # alignment
    if check_alg == 1 and len(set(map(len, result.values()))) != 1:
        sys.exit('\nERROR: Input sequence file must be an alignment.\n')
    return result

###---FUNCTION: Translate DNA sequence aligned at codon level---###
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

if __name__ == '__main__':

    # parse arguments
    params = argparse.ArgumentParser()
    params.add_argument('-i', '--input_file', required=True, type=str, help='Input evol-closest-pairs file to analyse')
    params.add_argument('-f', '--fa', required=True, type=str, help='FASTA alignment')
    #params.add_argument('--glyco_site', type=int, nargs='+', help='Association analysis of glycosylation status for site X (can be more than 1, absolute numbering)')
    params = params.parse_args()

    # parse fasta aln
    fdat = {k:translatedDNA(v) for k, v in parsefasta(params.fa).items()}

    # parse fo subtype
    try:
        subtype = re.search('(H3N2|pH1N1|H1N1pdm09|BVic|BYam)', params.input_file).group()
        aln_subtype = {'H3N2':'H3', 'pH1N1':'pH1abs', 'H1N1pdm09':'pH1abs', 'BVic':'B73', 'BYam':'B73'}[subtype]
    except:
        raise Exception('Can\'t parse for aln_subtype.')

    # copy structure file into current working directory
    # only H3 is chains = 2
    struc_file = {'H3N2':'4we8_Vic11_Repair.pdb', 'H1N1pdm09':'4jtv_Repair.pdb', 'BVic':'4fqm_Repair.pdb', 'BYam':'4m44_Repair.pdb'}[subtype]
    struc_dir = eu('~/Dropbox/age_final_submission/ageflu/files/')

    cmd = 'cp {}{} ./'.format(struc_dir, struc_file)
    subprocess.call(cmd, shell=True)

    # parse input file
    fhandle = filter(None, open(params.input_file, 'rU').readlines())

    all_mutation_list = []
    pt_to_charge_change_count = {}
    pt_to_mutation_count = {}
    pt_to_glyco_change_count = {}
    pos_to_pt_glycogainloss_count = {}

    for line in fhandle:
        line = line.strip().split('\t')
        pt = line[0]
        if not re.search('(AC|CA|AA|CC)', pt):
            continue

        anc, desc = line[1:3]
        try:
            mutations = line[3].split(',')
        except:
            mutations = []
        all_mutation_list = list(set(all_mutation_list + mutations))

        glyco_change_binary = 0

        if len(mutations) > 0 :
            for mutation in mutations:
                wt, mt = mutation[0], mutation[-1]

                # charge change
                charge_change_binary = 0
                if re.search('(R|H|K)', wt):
                    wt_charge = 1
                elif re.search('(D|E)', wt):
                    wt_charge = -1
                else:
                    wt_charge = 0

                if re.search('(R|H|K)', mt):
                    mt_charge = 1
                elif re.search('(D|E)', mt):
                    mt_charge = -1
                else:
                    mt_charge = 0

                charge_change = mt_charge - wt_charge
                if charge_change != 0:
                    charge_change_binary = 1

                # mutation count - stability
                try:
                    pt_to_mutation_count[pt][mutation] += 1
                except:
                    try:
                        pt_to_mutation_count[pt][mutation] = 1
                    except:
                        pt_to_mutation_count[pt] = {mutation:1}

                # glycosylation change
                #if glyco_change_binary == 0:
                abs_position_minus_1 = int(re.search('\d+', mutation).group()) - 1

                for abs_pos_iterator in xrange(abs_position_minus_1-2, abs_position_minus_1+1, 1):
                    anc_triplet = fdat[anc][abs_pos_iterator:abs_pos_iterator+3]
                    desc_triplet = fdat[desc][abs_pos_iterator:abs_pos_iterator+3]

                    # can't determine if there is a glycosylation change if there are any unknowns
                    if re.search('(X|-)', anc_triplet) or re.search('(X|-)', desc_triplet):
                        continue
                    else:
                        anc_glyco_binary = 1 if re.match('N[^P][ST]', anc_triplet) else 0
                        desc_glyco_binary = 1 if re.match('N[^P][ST]', desc_triplet) else 0

                        if anc_glyco_binary != desc_glyco_binary:
                            glyco_change_binary = 1

                            # gain glyco
                            if desc_glyco_binary == 1:
                                try:
                                    pos_to_pt_glycogainloss_count[abs_position_minus_1+1][pt]['1'] += 1
                                except:
                                    try:
                                        pos_to_pt_glycogainloss_count[abs_position_minus_1+1][pt]['1'] = 1
                                    except:
                                        try:
                                            pos_to_pt_glycogainloss_count[abs_position_minus_1+1][pt] = {'1':1}
                                        except:
                                            pos_to_pt_glycogainloss_count[abs_position_minus_1+1] = {pt:{'1':1}}
                            # lose glycosylation
                            else:
                                try:
                                    pos_to_pt_glycogainloss_count[abs_position_minus_1+1][pt]['0'] += 1
                                except:
                                    try:
                                        pos_to_pt_glycogainloss_count[abs_position_minus_1+1][pt]['0'] = 1
                                    except:
                                        try:
                                            pos_to_pt_glycogainloss_count[abs_position_minus_1+1][pt] = {'0':1}
                                        except:
                                            pos_to_pt_glycogainloss_count[abs_position_minus_1+1] = {pt:{'0':1}}
                            break

                # charge change
                if charge_change_binary > 0:
                    try:
                        pt_to_charge_change_count[pt]['1'] += 1
                    except:
                        try:
                            pt_to_charge_change_count[pt]['1'] = 1
                        except:
                            pt_to_charge_change_count[pt] = {'1': 1}
                else:
                    try:
                        pt_to_charge_change_count[pt]['0'] += 1
                    except:
                        try:
                            pt_to_charge_change_count[pt]['0'] = 1
                        except:
                            pt_to_charge_change_count[pt] = {'0':1}

            # glyco change
            if glyco_change_binary > 0:
                try:
                    pt_to_glyco_change_count[pt]['1'] += 1
                except:
                    try:
                        pt_to_glyco_change_count[pt]['1'] = 1
                    except:
                        pt_to_glyco_change_count[pt] = {'1':1}
            else:
                try:
                    pt_to_glyco_change_count[pt]['0'] += 1
                except:
                    try:
                        pt_to_glyco_change_count[pt]['0'] = 1
                    except:
                        pt_to_glyco_change_count[pt] = {'0':1}

    # stability change analyses

    # perform foldx using fold_stability_master_v2.py
    def perform_foldx(mut_list, p_index, q):
        # cp PDB files
        '''newpdb = '{}-proc{}.pdb'.format(re.sub('\.pdb$', '', struc_file), p_index)
        cmd = 'cp {} {}'.format(struc_file, newpdb)
        subprocess.call(cmd, shell=True)

        # perform foldx analyses
        mutfname = 'physico_mutation_list_proc{}.txt'.format(p_index)
        with open(mutfname, 'w') as output:
            for mut in mut_list:
                output.write('{_mut}\t{_mut}\n'.format(_mut=mut))'''

        foldx_outfname = 'foldx_output_proc{}.txt'.format(p_index)
        #if re.search('4we8', newpdb):
        '''cmd = ['python', eu('~/Dropbox/python_scripts/foldx_stability_master_v2.py'), '--mutfile', mutfname, '--pdb', newpdb, '--outfname', foldx_outfname, '--reference', aln_subtype, '--chains', '2']
        #else:
        #cmd = ['python', eu('~/Dropbox/python_scripts/foldx_stability_master_v2.py'), '--mutfile', mutfname, '--pdb', newpdb, '--outfname', foldx_outfname, '--reference', aln_subtype]
        foldx = subprocess.Popen(cmd)
        foldx.wait()'''

        fhandle = filter(None, [_.strip() for _ in open(foldx_outfname, 'rU')])
        fhandle.pop(0)

        foldxresults = {}
        for line in fhandle:
            line = line.strip().split('\t')
            main_sub_binary = line[1]

            if main_sub_binary == '*':
                mut = line[0]
                try:
                    ddG, ddG_sd = line[-2:]
                    foldxresults[mut] = round(float(ddG), 4)
                except:
                    continue # continue if no foldx result

        q.put(foldxresults)

    ncpu = mp.cpu_count()
    increment = int(round(len(all_mutation_list)/ncpu))
    Processes = []
    mut_to_foldxresults = {}

    # multi-proc setup
    manager = mp.Manager()

    # shared memory
    queue = manager.Queue()

    for p in xrange(ncpu):
        if p == ncpu - 1:
            curr_mut_list = all_mutation_list[p*increment:]
        else:
            curr_mut_list = all_mutation_list[p*increment:(p*increment) + increment]

        if len(curr_mut_list) == 0:
            continue

        proc = mp.Process(target=perform_foldx, args=(curr_mut_list, p, queue))
        proc.start()
        Processes.append(proc)
        time.sleep(5)

    # collect results to dictionary
    for p in xrange(len(Processes)):
        mut_to_foldxresults.update(queue.get())

    for proc in Processes:
        proc.join()

    # count
    pt_to_stability_change_count = {}

    for pt, mutation_to_count in pt_to_mutation_count.items():
        for mutation, count in mutation_to_count.items():
            if mutation in mut_to_foldxresults:
                if mut_to_foldxresults[mutation] > 0.46 or mut_to_foldxresults[mutation] < -0.46:
                    try:
                        pt_to_stability_change_count[pt]['1'] += count
                    except:
                        try:
                            pt_to_stability_change_count[pt]['1'] = count
                        except:
                            pt_to_stability_change_count[pt] = {'1':count}
                else:
                    try:
                        pt_to_stability_change_count[pt]['0'] += count
                    except:
                        try:
                            pt_to_stability_change_count[pt]['0'] = count
                        except:
                            pt_to_stability_change_count[pt] = {'0':count}

    # stats analyses and print output
    outfname = '{}_physicochemical_properties_association.txt'.format(subtype)

    with open(outfname, 'w') as output:
        output.write('Charge change\t\t\t\tStability change\t\t\t\tGlycosylation change\t\t\t\t\n'
                     'p-value\tOR\t90% CI\tn\tp-value\tOR\t90% CI\tn\tp-value\tOR\t90% CI\tn\n')

        # charge changes
        p11, p10, p01, p00 = get_contingency_counts(pt_to_charge_change_count)
        OR_uMLE, wald_MLE_CI90L, wald_MLE_CI90R, wald_MLE_CI95L, wald_MLE_CI95R, pval, ss_CI90L, ss_CI90R, ss_CI95L, ss_CI95R = calculateOR(p11, p10, p01, p00)
        n = sum([p11, p10, p01, p00])
        if n < 1000:
            try:
                CI90L = round(ss_CI90L, 2)
            except:
                CI90L = ss_CI90L

            try:
                CI90R = round(ss_CI90R, 2)
            except:
                CI90R = ss_CI95R
        else:
            try:
                CI90L = round(wald_MLE_CI90L, 2)
            except:
                CI90L = wald_MLE_CI90L
            try:
                CI90R = round(wald_MLE_CI90R, 2)
            except:
                CI90R = wald_MLE_CI90R

        try:
            pval = round(pval, 3)
        except:
            pval = pval

        try:
            OR_uMLE = round(OR_uMLE, 2)
        except:
            OR_uMLE = OR_uMLE

        output.write('{}\t{}\t({}, {})\t{}'.format(pval, OR_uMLE, CI90L, CI90R , n))


        # stability
        p11, p10, p01, p00 = get_contingency_counts(pt_to_stability_change_count)
        OR_uMLE, wald_MLE_CI90L, wald_MLE_CI90R, wald_MLE_CI95L, wald_MLE_CI95R, pval, ss_CI90L, ss_CI90R, ss_CI95L, ss_CI95R = calculateOR(p11, p10, p01, p00)
        n = sum([p11, p10, p01, p00])
        if n < 1000:
            try:
                CI90L = round(ss_CI90L, 2)
            except:
                CI90L = ss_CI90L

            try:
                CI90R = round(ss_CI90R, 2)
            except:
                CI90R = ss_CI95R
        else:
            try:
                CI90L = round(wald_MLE_CI90L, 2)
            except:
                CI90L = wald_MLE_CI90L
            try:
                CI90R = round(wald_MLE_CI90R, 2)
            except:
                CI90R = wald_MLE_CI90R

        try:
            pval = round(pval, 3)
        except:
            pval = pval

        try:
            OR_uMLE = round(OR_uMLE, 2)
        except:
            OR_uMLE = OR_uMLE

        output.write('\t{}\t{}\t({}, {})\t{}'.format(pval, OR_uMLE, CI90L, CI90R , n))

        # glycosylation changes
        p11, p10, p01, p00 = get_contingency_counts(pt_to_glyco_change_count)
        OR_uMLE, wald_MLE_CI90L, wald_MLE_CI90R, wald_MLE_CI95L, wald_MLE_CI95R, pval, ss_CI90L, ss_CI90R, ss_CI95L, ss_CI95R = calculateOR(p11, p10, p01, p00)
        n = sum([p11, p10, p01, p00])
        if n < 1000:
            try:
                CI90L = round(ss_CI90L, 2)
            except:
                CI90L = ss_CI90L

            try:
                CI90R = round(ss_CI90R, 2)
            except:
                CI90R = ss_CI95R
        else:
            try:
                CI90L = round(wald_MLE_CI90L, 2)
            except:
                CI90L = wald_MLE_CI90L
            try:
                CI90R = round(wald_MLE_CI90R, 2)
            except:
                CI90R = wald_MLE_CI90R

        try:
            pval = round(pval, 3)
        except:
            pval = pval

        try:
            OR_uMLE = round(OR_uMLE, 2)
        except:
            OR_uMLE = OR_uMLE

        output.write('\t{}\t{}\t({}, {})\t{}\n'.format(pval, OR_uMLE, CI90L, CI90R , n))

        # glycosylation gain/loss details for each position
        output.write('\npos\tCA\t\tAA\t\tAC\t\tCC\t\t\n')
        output.write('\tG\tL\tG\tL\tG\tL\tG\tL\n')
        pos_to_pt_glycogainloss_count = {}

        for pos in pos_to_pt_glycogainloss_count:
            try:
                CA_1 = pos_to_pt_glycogainloss_count[pos]['CA']['1']
            except:
                CA_1 = 0
            try:
                CA_0 = pos_to_pt_glycogainloss_count[pos]['CA']['0']
            except:
                CA_0 = 0
            try:
                AA_1 = pos_to_pt_glycogainloss_count[pos]['AA']['1']
            except:
                AA_1 = 0
            try:
                AA_0 = pos_to_pt_glycogainloss_count[pos]['AA']['0']
            except:
                AA_0 = 0
            try:
                AC_1 = pos_to_pt_glycogainloss_count[pos]['AC']['1']
            except:
                AC_1 = 0
            try:
                AC_0 = pos_to_pt_glycogainloss_count[pos]['AC']['0']
            except:
                AC_0 = 0
            try:
                CC_1 = pos_to_pt_glycogainloss_count[pos]['CC']['1']
            except:
                CC_1 = 0
            try:
                CC_0 = pos_to_pt_glycogainloss_count[pos]['CC']['0']
            except:
                CC_0 = 0
            output.write('{}\n'.format('\t'.join(map(str, [pos, CA_1, CA_0, AA_1, AA_0, AC_1, AC_0, CC_1, CC_0]))))

        """
        # glycosylation status association analysis
        if params.glyco_site:
            output.write('Pos\tchild_glyco0\tchild_glyco1\tadult_glyco0\tadult_glyco1\tOR\twald90CI\tss90CI\tpval\n')
            for glycopos in params.glyco_site:

                agetype_to_glycostatus_count = {}

                for header, sequence in fdat.items():
                    glyco_triplet = sequence[glycopos-1:glycopos+2]
                    # continue if there are any unknown residue
                    if re.search('(X|-)', glyco_triplet):
                        continue

                    glyco_binary = 1 if re.match('N[^P][ST]', glyco_triplet) else 0

                    # classify based on age threshold
                    age = int(re.search('AGE(\d+)', header).group(1))
                    if 0 <= age <= params.age[0]:
                        agetype = 'C'
                    elif params.age[-1] <= age <= 120:
                        agetype = 'A'
                    else:
                        continue

                    # store count
                    try:
                        agetype_to_glycostatus_count[agetype][glyco_binary] += 1
                    except:
                        try:
                            agetype_to_glycostatus_count[agetype][glyco_binary] = 1
                        except:
                            agetype_to_glycostatus_count[agetype] = {glyco_binary:1}

                try:
                    C_0 = agetype_to_glycostatus_count['C'][0]
                except:
                    C_0 = 0

                try:
                    C_1 = agetype_to_glycostatus_count['C'][1]
                except:
                    C_1 = 0

                try:
                    A_0 = agetype_to_glycostatus_count['A'][0]
                except:
                    A_0 = 0

                try:
                    A_1 = agetype_to_glycostatus_count['A'][1]
                except:
                    A_1 = 0

                OR_uMLE, wald_MLE_CI90L, wald_MLE_CI90R, wald_MLE_CI95L, wald_MLE_CI95R, pval, ss_CI90L, ss_CI90R, ss_CI95L, ss_CI95R = calculateOR(C_0, C_1, A_0, A_1)
                output.write('{}\n'.format('\t'.join(map(str, [glycopos, C_0, C_1, A_0, A_1, OR_uMLE, '{},{}'.format(wald_MLE_CI90L, wald_MLE_CI90R), '{},{}'.format(ss_CI90L, ss_CI90R), pval]))))"""


    exit(0)