#############################################
#  Python script to parse closest pairs     #
#############################################

# parse fasta file
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

# translate codon-aligned nucleotide sequences
def translatedDNA(dna):
    protein = []
    for _ in xrange(0,len(dna),3):
        codon = dna[_:_+3]
        if codon in dnacodontable:
            if dnacodontable[codon] == 'stop': break
            protein.append(dnacodontable[codon])
        elif re.match('---',codon): protein.append('-')
        else: protein.append('X')
    return ''.join(protein)

# get pairwise substitutions
def get_pairwise_substitutions(anc_seq, desc_seq):
    # nucleotide sequences input
    if set(list(anc_seq)) <= set(list('atgcn-')):
        transitions, transversions, nonsyn_sub, syn_sub = 0, 0, 0, 0
        for _ in xrange(0, len(anc_seq), 3):
            anc_codon = anc_seq[_:_+3]
            desc_codon = desc_seq[_:_+3]
            if anc_codon == desc_codon:
                continue
            else:
                unknown_res_binary = 0
                nuc_sub = 0
                for j in xrange(3):
                    if re.match('(n|-)', anc_codon[j]) or re.match('(n|-)', desc_codon[j]):
                        unknown_res_binary = 1
                        continue
                    elif anc_codon[j] != desc_codon[j]:
                        if re.search('(ag|ga|ct|tc)', ''.join([anc_codon[j], desc_codon[j]])):
                            transitions += 1
                            nuc_sub += 1
                        else:
                            transversions += 1
                            nuc_sub += 1

                if unknown_res_binary == 1:
                    continue
                if dnacodontable[anc_codon] != dnacodontable[desc_codon]:
                    nonsyn_sub += nuc_sub
                else:
                    syn_sub += nuc_sub
        return {'transitions':transitions, 'transversions':transversions, 'nonsyn':nonsyn_sub, 'syn':syn_sub}
    # protein sequences input
    else:
        mutation_list = []
        for _ in xrange(len(anc_seq)):
            if re.search('(-|X)', anc_seq[_]) or re.search('(-|X)', desc_seq[_]):
                continue
            elif anc_seq[_] != desc_seq[_]:
                mutation_list.append('{}{}{}'.format(anc_seq[_], _+1, desc_seq[_]))
        return mutation_list

def get_least_unknown_residue_sequences(seq_dict, dist_dict_to_edit={}):
    seqid_to_unknownrescount = {seqid:len(re.findall('(-|n|X)', sequence, re.I)) for seqid, sequence in seq_dict.items()}
    seqid_with_least_unknown_res = [seqid for seqid, count in seqid_to_unknownrescount.items() if count == min(seqid_to_unknownrescount.values())]
    if len(dist_dict_to_edit) > 0:
        for seqid in list(set(seq_dict.keys())-set(seqid_with_least_unknown_res)):
            del dist_dict_to_edit[seqid]
        return seqid_with_least_unknown_res, dist_dict_to_edit
    else:
        return seqid_with_least_unknown_res

if __name__ == '__main__':
    # global parameters
    dnacodontable = {'ttt':'F', 'ttc':'F', 'tta':'L', 'ttg':'L', 'ctt':'L', 'ctc':'L', 'cta':'L', 'ctg':'L', 'att':'I', 'atc':'I', 'ata':'I', 'atg':'M', 'gtt':'V', 'gtc':'V', 'gta':'V', 'gtg':'V', 'tct':'S', 'tcc':'S', 'tca':'S', 'tcg':'S', 'cct':'P', 'ccc':'P', 'cca':'P', 'ccg':'P', 'act':'T', 'acc':'T', 'aca':'T', 'acg':'T', 'gct':'A', 'gcc':'A', 'gca':'A', 'gcg':'A', 'tat':'Y', 'tac':'Y', 'taa':'stop', 'tag':'stop', 'cat':'H', 'cac':'H', 'caa':'Q', 'cag':'Q', 'aat':'N', 'aac':'N', 'aaa':'K', 'aag':'K', 'gat':'D', 'gac':'D', 'gaa':'E', 'gag':'E', 'tgt':'C', 'tgc':'C', 'tga':'stop', 'tgg':'W', 'cgt':'R', 'cgc':'R', 'cga':'R', 'cgg':'R', 'agt':'S', 'agc':'S', 'aga':'R', 'agg':'R', 'ggt':'G', 'ggc':'G', 'gga':'G', 'ggg':'G'}

    # parameters to obtain reference positions
    from os.path import expanduser
    import re, sys, itertools
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

    # parse argument
    import argparse
    params = argparse.ArgumentParser(description='Parse tree for evolutionarily closest pairs.')
    params.add_argument('-t', '--tree', type=str, required=True, help='rooted NEWICK tree file.')
    params.add_argument('-f', '--aln', type=str, required=True, help='Codon-aligned nucleotide fasta of sequences.')
    params.add_argument('-n', '--nr', type=str, required=True, help='Non-redundancy CD-HIT clstr file.')
    params.add_argument('-c', '--child', type=int, nargs=2 , default=[0, 5], help='Children age range (default: %(default)s).')
    params.add_argument('-a', '--adult', type=int, nargs=2 , default=[35, 120], help='Adult age range (default: %(default)s).')
    params.add_argument('-mm', '--maxmut', type=int, default=5, help='Maximum number of amino acid substitutions of pair (default: %(default)s).')
    params.add_argument('-pd', '--patdist', type=float, default=0.007, help='Maximum patristic distance of pair (default: %(default)s).')
    params.add_argument('-o', '--order', type=int, default=3, help='Analyze up to the n-th ordinal set of pairs (default: %(default)s).')
    params = params.parse_args()

    # query subtype analysed
    try:
        query_subtype = re.search('(H3N2|pH1N1|H1N1pdm09|BVic|BYam)', params.tree).group().upper()
        if query_subtype == 'PH1N1':
            query_subtype = 'H1N1PDM09'
    except:
        query_subtype = raw_input('\nWARNING: Can\'t parsed query subtype from tree name. Enter subtype (H3N2/pH1N1/BVic/BYam): ')

    # output filename
    outfname = '{}C{}_{}A{}_PD{}_MM{}_CL{}_{}'.format(params.child[0], params.child[-1], params.adult[0], params.adult[-1], params.patdist, params.maxmut, params.order, re.sub('[^/]*/', '', params.tree))

    # parse fasta alignment
    print ('\n...parsing fasta aln...')
    fdat_nuc = parsefasta(params.aln)
    fdat_aa = {k:translatedDNA(v) for k,v in fdat_nuc.items()}
    isolateid_to_fdatheader = {re.search('(GBISL|EPIISL|GB_ISL_|EPI_ISL_)\d+', header).group().replace('_', ''):header for header in fdat_nuc.keys()}
    fdatheader_to_isolateid = {v:k for k,v in isolateid_to_fdatheader.items()}
    isolateid_to_age = {isolateid:int(re.search('AGE(\d+)', header).group(1)) for isolateid, header in isolateid_to_fdatheader.items()}

    with open('dates.txt', 'w') as output:
        for isolateid, header in isolateid_to_fdatheader.items():
            date = re.search('_(\d+\.\d+)_', header).group(1)
            output.write('{},{}\n'.format(isolateid, date))
    exit(1)
    # sequences to ignore (outside of child_min and adult_max age)
    isolates_to_ignore = [isolateid for isolateid, age in isolateid_to_age.items() if age < params.child[0] or age > params.adult[-1]]

    # parse clstr file and sort each identical sequence cluster (>1 member sequence) into children/adult age categories
    print ('...parsing cd-hit clstr file...')
    isolateid_to_representatives, isolateid_to_agecategory = {}, {} # isolates_to_ignore include those of ambiguous age categories + 0 pat dist sequences with > minimal number of n/-
    fhandle = filter(None, open(params.nr, 'rU').readlines())
    for key, group in itertools.groupby(fhandle, lambda _: re.search('^>Cluster', _)):
        if not key:
            cluster = list(group)
            if len(cluster) > 1:
                cluster = [re.search('(GBISL|EPIISL|GB_ISL_|EPI_ISL_)\d+', header).group().replace('_', '') for header in cluster]
                cluster_age = [isolateid_to_age[isolateid] for isolateid in cluster]
                # all children or adults sequence clusters
                if all([params.child[0] <= age <= params.child[-1] for age in cluster_age]):
                    isolateid_to_agecategory[isolateid] = 'C'
                    for isolateid in cluster:
                        isolateid_to_representatives[isolateid] = cluster
                elif all([params.adult[0] <= age <= params.adult[-1] for age in cluster_age]):
                    isolateid_to_agecategory[isolateid] = 'A'
                    for isolateid in cluster:
                        isolateid_to_representatives[isolateid] = cluster
                # uncategorized age present in sequence cluster
                elif all([params.child[0] <= age <= params.adult[-1] for age in cluster_age]) and any([params.child[-1] < age < params.adult[0] for age in cluster_age]):
                    isolateid_to_agecategory[isolateid] = 'U'
                    for isolateid in cluster:
                        isolateid_to_representatives[isolateid] = [member for _, member in enumerate(cluster) if params.child[-1] < cluster_age[_] < params.adult[0]]
                else:
                    isolates_to_ignore = list(set(isolates_to_ignore)|set(cluster))

    # parse tree
    print ('...parsing tree...')
    import ete3
    try:
        tree = ete3.Tree(params.tree)
    except:
        sys.exit('\nERROR: Unable to parse tree using ete3.\n')
    tree.ladderize()
    # get nodes to isolate ids
    leaf_to_isolateid = {leaf:re.search('(EPI_ISL_|EPIISL|GB_ISL_|GBISL)\d+', leaf.get_leaf_names()[0]).group().replace('_', '') for leaf in tree.get_leaves()}
    isolateid_to_leaf = {v:k for k,v in leaf_to_isolateid.items()}

    import random
    from decimal import *
    print ('...finding closest pairs...')
    anc_to_desc = {}
    desc_to_anc = {}
    pair_to_subinfo = {}
    pair_to_pwdist = {}
    # traverse through tree level-order
    for node in tree.traverse():
        if node.is_leaf():
            continue

        # get children leaf nodes and their distances to current node
        childleaf_to_distance = {child:Decimal(child.get_distance(node)).quantize(Decimal('1e-4')) for child in node.children if child.is_leaf()}

        if len(childleaf_to_distance) > 0:
            # closest child leaf from node - decide closest child leaf using patristic distances without truncating at 1e-5
            nearest_childleaf = [child for child in childleaf_to_distance.keys() if childleaf_to_distance[child] == min(childleaf_to_distance.values())]
            # remove sequences in cases where distance = 0 and aren't the least number of unknown residues
            if min(childleaf_to_distance.values()) == 0 and len(nearest_childleaf) > 1:
                nearest_childleaf, childleaf_to_distance = get_least_unknown_residue_sequences({child:fdat_nuc[isolateid_to_fdatheader[leaf_to_isolateid[child]]] for child in nearest_childleaf}, dist_dict_to_edit=childleaf_to_distance)

            # define nearest child leaf/leaves
            nearest_childleaf_distance_to_node = {child:childleaf_to_distance[child] for child in nearest_childleaf}
            for child in nearest_childleaf:
                del childleaf_to_distance[child]

            # get grandchildren leaves
            descendant_internal_nodes = [child for child in node.get_descendants() if not child.is_leaf() and node.get_distance(child) < params.patdist]
            if len(descendant_internal_nodes) > 0:
                for childnode in descendant_internal_nodes:
                    grandchildren_leaves = [grandchild for grandchild in childnode.children if grandchild.is_leaf()]
                    if len(grandchildren_leaves) > 0:
                        grandchildleaf_to_distance = {grandchild:Decimal(grandchild.get_distance(childnode)).quantize(Decimal('1e-4')) for grandchild in grandchildren_leaves}
                        nearest_grandchild_with_zero_distance = [grandchild for grandchild in grandchildren_leaves if grandchildleaf_to_distance[grandchild] == 0]
                        # remove sequences in cases where distance = 0 and aren't the least number of unknown residues
                        # we do this because we don't want to re-count possible identical strains
                        if len(nearest_grandchild_with_zero_distance) > 1:
                            nearest_grandchild_with_zero_distance_to_keep = get_least_unknown_residue_sequences({grandchild:fdat_nuc[isolateid_to_fdatheader[leaf_to_isolateid[grandchild]]] for grandchild in nearest_grandchild_with_zero_distance})
                            for grandchild in list(set(nearest_grandchild_with_zero_distance)-set(nearest_grandchild_with_zero_distance_to_keep)):
                                grandchildren_leaves.remove(grandchild)
                        # update grandchild distance to node
                        childleaf_to_distance.update({grandchild:Decimal(grandchild.get_distance(node)).quantize(Decimal('1e-4')) for grandchild in grandchildren_leaves})

            # if there are no children/grandchildren left for pairing
            if len(childleaf_to_distance) == 0:
                continue

            for anc_leaf in nearest_childleaf:
                # ignore ambiguous age isolates
                anc = leaf_to_isolateid[anc_leaf]
                if anc in isolates_to_ignore:
                    continue

                # find closest pairing descendants
                # get pairwise distance to ancestor
                childleaf_to_pwdistance = {}
                for child, distance in childleaf_to_distance.items():
                    pwdistance = nearest_childleaf_distance_to_node[anc_leaf] + distance
                    childleaf_to_pwdistance[child] = pwdistance
                    pair_to_pwdist[(anc, leaf_to_isolateid[child])] = pwdistance

                # analyze up to the n-th order set of pairs
                for pwdistance in sorted(set(childleaf_to_pwdistance.values()))[:params.order]:
                    # patristic distance threshold
                    if pwdistance > params.patdist:
                        continue

                    children_with_distance = [child for child, child_distance in childleaf_to_pwdistance.items() if child_distance == pwdistance]
                    for desc_leaf in children_with_distance:

                        desc = leaf_to_isolateid[desc_leaf]
                        if desc in isolates_to_ignore:
                            continue

                        # descendant should be further away from node
                        if childleaf_to_distance[desc_leaf] <= nearest_childleaf_distance_to_node[anc_leaf]:
                            continue

                        substitutions = get_pairwise_substitutions(fdat_aa[isolateid_to_fdatheader[anc]], fdat_aa[isolateid_to_fdatheader[desc]])
                        # amino acid substitutions threshold
                        if len(substitutions) > params.maxmut:
                            continue

                        # check that there is no nearer ancestor to descendant
                        if desc in desc_to_anc:
                            prev_anc_list = desc_to_anc[desc]
                            prev_anc_to_pwdist = {prev_anc:pair_to_pwdist[(prev_anc, desc)] for prev_anc in prev_anc_list}
                            prev_anc_to_pwdist[anc] = pwdistance
                            minimum_distance = min(prev_anc_to_pwdist.values())

                            for prev_anc, prev_anc_pwdist in prev_anc_to_pwdist.items():
                                if prev_anc_pwdist == minimum_distance:
                                    if prev_anc in anc_to_desc:
                                        if desc not in anc_to_desc[prev_anc]:
                                            anc_to_desc[prev_anc].append(desc)
                                    else:
                                        anc_to_desc[prev_anc] = [desc]

                                    if desc in desc_to_anc:
                                        if prev_anc not in desc_to_anc[desc]:
                                            desc_to_anc[desc].append(prev_anc)
                                    else:
                                        desc_to_anc[desc] = [prev_anc]

                                    if (prev_anc, desc) not in pair_to_subinfo:
                                        pair_to_subinfo[(prev_anc, desc)] = (substitutions, get_pairwise_substitutions(fdat_nuc[isolateid_to_fdatheader[anc]], fdat_nuc[isolateid_to_fdatheader[desc]]))
                                else:
                                    try:
                                        anc_to_desc[prev_anc].remove(desc)
                                        if len(anc_to_desc[prev_anc]) == 0:
                                            del anc_to_desc[prev_anc]
                                    except:
                                        pass

                                    try:
                                        desc_to_anc[desc].remove(prev_anc)
                                        if len(desc_to_anc[desc]) == 0:
                                            del desc_to_anc[desc]
                                    except:
                                        pass
                        else:
                            try:
                                anc_to_desc[anc].append(desc)
                            except:
                                anc_to_desc[anc] = [desc]

                            try:
                                desc_to_anc[desc].append(anc)
                            except:
                                desc_to_anc[desc] = [anc]

                            pair_to_subinfo[(anc, desc)] = (substitutions, get_pairwise_substitutions(fdat_nuc[isolateid_to_fdatheader[anc]], fdat_nuc[isolateid_to_fdatheader[desc]]))

    # write closest pairs output
    print ('...writing output...')
    with open('ageflu_evol-closest-pairs.{}.txt'.format(outfname), 'w') as output:
        output.write('pair_type\tanc\tdesc\tmutation(ABS)\tmutation(REF)\n')
        for anc, desclist in anc_to_desc.items():
            for desc in desclist:

                # determine pair type
                if params.child[0] <= isolateid_to_age[anc] <= params.child[-1]:
                    anc_age_type = 'C'
                elif params.adult[0] <= isolateid_to_age[anc] <= params.adult[-1]:
                    anc_age_type = 'A'
                else:
                    anc_age_type = 'N'

                if params.child[0] <= isolateid_to_age[desc] <= params.child[-1]:
                    desc_age_type = 'C'
                elif params.adult[0] <= isolateid_to_age[desc] <= params.adult[-1]:
                    desc_age_type = 'A'
                else:
                    desc_age_type = 'N'

                if anc_age_type == 'N' or desc_age_type == 'N':
                    pair_age_type = 'NC'
                else:
                    pair_age_type = '{}{}'.format(anc_age_type, desc_age_type)

                # get substitution info
                substitutions, nuc_sub_info = pair_to_subinfo[(anc, desc)]
                substitutions_abs_positions = [int(re.search('\d+', _).group()) for _ in substitutions]

                # get reference positions
                substitutions_ref_positions = [queryST_to_refST_to_AbNum_to_RefNum[query_subtype][ref_subtype_dictionary[query_subtype]][abs_pos] if abs_pos in queryST_to_refST_to_AbNum_to_RefNum[query_subtype][ref_subtype_dictionary[query_subtype]] else '-' for abs_pos in substitutions_abs_positions]

                # write to output
                write_line = map(str, [pair_age_type, isolateid_to_fdatheader[anc], isolateid_to_fdatheader[desc], ','.join(substitutions), ','.join(map(str, substitutions_ref_positions))])
                output.write('{}\n'.format('\t'.join(write_line)))
                output.flush()

    print ('...done.\n')
    exit(0)
