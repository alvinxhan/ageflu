#!/usr/bin/env python

###---FUNCTION: parse pdb file for protein structure coordinates---###
def parsepdb(file):
    fhandle = filter(None, open(file, 'rU').readlines())

    # amino acid
    aadict = {'ARG':'R', 'HIS':'H', 'LYS':'K', 'ASP':'D', 'GLU':'E', 'SER':'S', 'THR':'T', 'ASN':'N', 'GLN':'Q', 'CYS':'C',
              'GLY':'G', 'PRO':'P', 'ALA':'A', 'VAL':'V', 'ILE':'I', 'LEU':'L', 'MET':'M', 'PHE':'F', 'TYR':'Y', 'TRP':'W'}

    chain_to_seqindex_to_pos_res = {}
    chain_to_seqindex = {}
    for line in fhandle:
        # all info will be based on alpha C
        if re.search("^ATOM", line) and line[12:16].replace(" ", "") == "CA":
            structural_pos = int(line[22:26])
            res = line[17:20]
            chain = line[21]

            try:
                seqindex = chain_to_seqindex[chain]
            except:
                seqindex = 0
                chain_to_seqindex[chain] = seqindex

            try:
                chain_to_seqindex_to_pos_res[chain][seqindex] = (structural_pos, aadict[res])
            except:
                chain_to_seqindex_to_pos_res[chain] = {seqindex:(structural_pos, aadict[res])}

            chain_to_seqindex[chain] += 1

    chain_to_sequence = {}
    for chain in chain_to_seqindex_to_pos_res.keys():
        chain_to_sequence[chain] = ''.join([chain_to_seqindex_to_pos_res[chain][seqindex][-1] for seqindex in chain_to_seqindex_to_pos_res[chain]])

    return chain_to_seqindex_to_pos_res, chain_to_sequence

###---FUNCTION: PARSE FASTA---###
def parsefasta(file, record_header_order=0):
    # check if file is of FASTA format
    fhandle = filter(None, open(file, "rU").readlines())
    if len(fhandle) == 0:
        print file
    if not re.search("^>", fhandle[0]):
        sys.exit()

    header_to_sequence = {}
    header_order = []
    sequence_len = 0
    for key, group in itertools.groupby(fhandle, lambda _: re.search("^>",_)):
        if key:
            header = group.next().strip().replace(">","")
        else:
            sequence = "".join(map(lambda _:_.strip(),list(group))).upper()
            # must be alignment
            if sequence_len == 0:
                sequence_len = len(sequence)
            elif len(sequence) != sequence_len:
                print ('\nERROR: FASTA is not an alignment.\n')
                sys.exit()
            header_to_sequence[header] = sequence
            if record_header_order == 1:
                header_order.append(header)

    if record_header_order == 1:
        return header_to_sequence, header_order
    else:
        return header_to_sequence

def get_ref_to_structural_pos(query_chain, fdat):

    result_dict = {}
    refseq = fdat['reference']
    structuralseq = fdat[query_chain]
    for i, refaa in enumerate(refseq):
        if refaa == '-' or structuralseq[i] == '-':
            continue
        refpos = len(refseq[:i+1].replace('-', ''))
        structuralseq_index = len(structuralseq[:i+1].replace('-', ''))-1
        try:
            structuralpos, structuralres = chain_to_seqindex_to_pos_res[query_chain][structuralseq_index]
        except:
            continue

        if structuralres != structuralseq[i]:
            print refpos, structuralpos, structuralres, structuralseq[i]
            print refseq, '\n', structuralseq, i
            sys.exit('\nERROR in alignment.\n')

        result_dict[refpos] = (structuralpos, structuralres)

    return result_dict

# perform foldx
def foldx(input_list, no_of_runs, pdbfname_to_use, ion, pH, temp, output_binary):
    indv_list_fname = 'individual_list_{}.txt'.format(pdbfname_to_use)
    with open(indv_list_fname, 'w') as output:
        if isinstance(input_list, list):
            output.write('{};\n'.format(','.join(input_list)))
        else:
            output.write('{};\n'.format(input_list))

    cmd = [expanduser('~/Dropbox/softwares/foldx_{}/foldx'.format(system)), '--command=BuildModel',
           '--pdb={}.pdb'.format(pdbfname_to_use), '--mutant-file={}'.format(indv_list_fname), '--ionStrength={}'.format(ion),
           '--pH={}'.format(pH), '--water=CRYSTAL', '--vdwDesign=2', '--out-pdb=true', '--pdbHydrogens=false',
           '--numberOfRuns={}'.format(no_of_runs), '--temperature={}'.format(temp),
           '--rotabaseLocation={}/Dropbox/softwares/foldx_{}/rotabase.txt'.format(home, system)]

    subprocess.call(cmd, stderr=subprocess.STDOUT, stdout=stdout_file) #stderr=subprocess.PIPE, stdout=subprocess.PIPE)

    if output_binary == 1:
        result_fhandle = filter(None, open('Raw_{}.fxout'.format(pdbfname_to_use), 'rU').readlines())
        record_result_BIN = 0
        run_to_result = {}
        for line in result_fhandle:
            if re.search('^Pdb\ttotal energy',line):
                record_result_BIN = 1
                continue
            if record_result_BIN > 0:
                if no_of_runs > 1:
                    run_index = int(re.search('(\d+)\.pdb', line).group(1)) + 1
                else:
                    run_index = int(re.search('(\d+)\.pdb', line).group(1))

                # get total energy
                energy = float(line.split('\t')[1])
                # check if this is WT or MT energy
                if re.search('^WT_', line):
                    res_type = 'WT'
                else:
                    res_type = 'MT'
                # save result
                try:
                    run_to_result[run_index][res_type] = energy
                except:
                    run_to_result[run_index] = {res_type:energy}

    # remove foldx output files and individual_list.txt
    #subprocess.call('rm WT_{}*.pdb individual_list.txt'.format(pdbfname_to_use), shell=True, stderr=subprocess.PIPE, stdout=subprocess.PIPE)
    subprocess.call('rm WT_{_pdbfname}*.pdb *{_pdbfname}*.fxout {_indv_list_fname}'.format(_pdbfname=pdbfname_to_use, _indv_list_fname=indv_list_fname), shell=True, stderr=subprocess.STDOUT, stdout=stdout_file)

    if output_binary == 1:
        return run_to_result
    else:
        return

def merge_two_dicts(x, y):
    z = x.copy()   # start with x's keys and values
    z.update(y)    # modifies z with y's keys and values & returns None
    return z

if __name__ == '__main__':
    # preset reference sequences
    #preset_reference_sequences = {'pH1abs' : 'MKAILVVLLYTFATANADTLCIGYHANNSTDTVDTVLEKNVTVTHSVNLLEDKHNGKLCKLRGVAPLHLGKCNIAGWILGNPECESLSTASSWSYIVETSSSDNGTCYPGDFIDYEELREQLSSVSSFERFEIFPKTSSWPNHDSNKGVTAACPHAGAKSFYKNLIWLVKKGNSYPKLSKSYINDKGKEVLVLWGIHHPSTSADQQSLYQNADAYVFVGSSRYSKKFKPEIAIRPKVRDQEGRMNYYWTLVEPGDKITFEATGNLVVPRYAFAMERNAGSGIIISDTPVHDCNTTCQTPKGAINTSLPFQNIHPITIGKCPKYVKSTKLRLATGLRNVPSIQSRGLFGAIAGFIEGGWTGMVDGWYGYHHQNEQGSGYAADLKSTQNAIDEITNKVNSVIEKMNTQFTAVGKEFNHLEKRIENLNKKVDDGFLDIWTYNAELLVLLENERTLDYHDSNVKNLYEKVRSQLKNNAKEIGNGCFEFYHKCDNTCMESVKNGTYDYPKYSEEAKLNREEIDGVKLESTRIYQILAIYSTVASSLVLVVSLGAISFWMCSNGSLQCRICI', 'H3' : 'QDLPGNDNSTATLCLGHHAVPNGTLVKTITDDQIEVTNATELVQSSSTGKICNNPHRILDGIDCTLIDALLGDPHCDVFQNETWDLFVERSKAFSNCYPYDVPDYASLRSLVASSGTLEFITEGFTWTGVTQNGGSNACKRGPGSGFFSRLNWLTKSGSTYPVLNVTMPNNDNFDKLYIWGIHHPSTNQEQTSLYVQASGRVTVSTRRSQQTIIPNIGSRPWVRGLSSRISIYWTIVKPGDVLVINSNGNLIAPRGYFKMRTGKSSIMRSDAPIDTCISECITPNGSIPNDKPFQNVNKITYGACPKYVKQNTLKLATGMRNVPEKQTRGLFGAIAGFIENGWEGMIDGWYGFRHQNSEGTGQAADLKSTQAAIDQINGKLNRVIEKTNEKFHQIEKEFSEVEGRIQDLEKYVEDTKIDLWSYNAELLVALENQHTIDLTDSEMNKLFEKTRRQLRENAEEMGNGCFKIYHKCDNACIESIRNGTYDHDVYRDEALNNRFQIKGVELKSGYKDWILWISFAISCFLLCVVLLGFIMWACQRGNIRCNICI','PB2' : 'MERIKELRDLMSQSRTREILTKTTVDHMAIIKKYTSGRQEKNPALRMKWMMAMKYPITADKRIIEMIPERNEQGQTLWSKTNDAGSDRVMVSPLAVTWWNRNGPTTSAVHYPKVYKTYFEKVERLKHGTFGPVHFRNQVKIXRRVDINPGHADLSAKEAQDVIMEVVFPNEVGARILTSESQLTITKEKKEELQDCKIAPLMVAYMLERELVRKTRFLPVAGGTSSVYIEVLHLTQGTCWEQMYTPGGEVRNDDVDQSLIIAARNIVRRATVSADPLASLLEMCHSTQIGGTRMVDILRQNPTEEQAVDICKAAMGLRISSSLSFGGFTFKRTSGSSVTKEEEVLTGNLQTLKIRVHEGYEEFTMVGRRATAILRKATRRLIQLIVSGRDEQSIAEAIIVAMVFSQEDCMIKAVRGDLNFVNRANQRLNPMHQLLRHFQKDAKVLFQNWGIEPIDNVMGMIGILPDMTPSTEMSLRGVRVSKMGVDEYSSTERVVVSIDRFLRVRDQRGNVLLSPEEVSETQGIERLTITYSSSMMWEINGPESVLVNTYQWIIRNWETVKIQWSQDPTMLYNKMEFEPFQSLVPKAARGQYSGFVRTLFQQMRDVLGTFDTVQIIKLLPFAAAPPEQSRMQFSSLTVNVRGSGMRILVRGNSPVFNYNKATKRLTVLGKDAGALTEDPDEGTAGVESAVLRGFLILGKEDKRYGPALSINELSNLAKGEKANVLIGQGDVVLVMKRKRDSSILTDSQTATKRIRMAIN', 'PB1' : 'MDVNPTLLFLKVPVQNAISTTFPYTGDPPYSHGTGTGYTMDTVNRTHQYSEKGKWTTNTETGAPQLNPIDGPLPEDNEPSGYAQTDCVLEAMAFLEESHPGIFENSCLETMEIVQQTRVDKLTQGRQTYDWTLNRNQPAATALANTIEIFRSNGLTANESGRLIDFLKDVMESMDKEEMEITTHFQRKRRVRDNMTKKMVTQRTIGKKKQRLNKRSYLIRALTLNTMTKDAERGKLKRRAIATPGMQIRGFVYFVETLARSICEKLEQSGLPVGGNEKKAKLANVVRKMMTNSQDTELSFTITGDNTKWNENQNPRMFLAMITYITRNQPEWFRNVLSIAPIMFSNKMARLGKGYMFESKSMKLRTQIPAEMLASIDLKYFNELTKKKIEKIRPLLIDGAASLSPGMMMGMFNMLSTVLGVSILNLGQKRYTKTTYWWDGLQSSDDFALIVNAPNHEGIQAGVDRFYRTCKLVGINMSKKKSYINRTGTFEFTSFFYRYGFVANFSMELPSFGVSGINESADMSIGVTVIKNNMINNDLGPATAQMALQLFIKDYRYTYRCHRGDTQIQTRRSFELKKLWEQTRSKAGLLVSDGGPNLYNIRNLHIPEVCLKWELMDEDYQGRLCNPLNPFVSHKEIESVNNAVVMPAHGPAKSMEYDAVATTHSWIPKRNRSILNTSQRGILEDEQMYQKCCNLFEKFFPSSSYRRPVGISSMVEAMVSRARIDARIDFESGRIKKEEFAEIMKICSTIEELRRQK', 'PA' : 'MEDFVRQCFNPMIVELAEKAMKEYGEDPKIETNKFAAICTHLEVCFMYSDFHFIDERSESIIVESGDPNALLKHRFEIIEGRDRTMAWTVVNSLCNTTGVEKPKFLPDLYDYKENRFIEIGVTRREVHTYYLEKANKIKSEKTHIHIFSFTGEEMATKADYTIDEESRARIKTRLFTIRQEMASRGLWDSFRQSERGEETIEEKFEITGTMRRLADQSLPPNFSSLENFRAYVDGFEPNGCIEGKLSQMSKEVNARIEPFLKTTPRPLRLPDGPPCSQRSKFLLMDALKLSIEDPSHEGEGIPLYDAIKCMKTFFGWKEPNIVKPHEKGINPNYLLAWKQVLAELQDIENEEKIPKTKNMKKTSQLKWALGENMAPEKVDFEDCKDVSDLTQYNSDEPESRSLASWIQSEFNKACELTDSSWIELDEIGEDVAPIEHIASMRRNYFTAEVSHCRATEYIMKGVYINTALLNASCAAMDDFQLIPMISKCRTKEGRRKTNLYGFIIKGRSHLRNDTDVVNYVSMEFSLTDPRLEPHKWEKYCVLEIGDMLLRTAVGQVSRPMFLYVRTNGTSKIKMKWGMEMRRCLLQSLQQIESMIEAESSVKEKDMTKEFFENKSETWPIGESPKGVEEGSIGKVCRTLLAKSVFNSLYASSQLEGFSAESRKLLLIAQALRDNLEPGTFDLGGLYEAIEECLINDPWVLLNASWFNSFLAHALK'}
    preset_reference_sequences = {'pH1abs' : 'MKAILVVLLYTFATANADTLCIGYHANNSTDTVDTVLEKNVTVTHSVNLLEDKHNGKLCKLRGVAPLHLGKCNIAGWILGNPECESLSTASSWSYIVETSSSDNGTCYPGDFIDYEELREQLSSVSSFERFEIFPKTSSWPNHDSNKGVTAACPHAGAKSFYKNLIWLVKKGNSYPKLSKSYINDKGKEVLVLWGIHHPSTSADQQSLYQNADAYVFVGSSRYSKKFKPEIAIRPKVRDQEGRMNYYWTLVEPGDKITFEATGNLVVPRYAFAMERNAGSGIIISDTPVHDCNTTCQTPKGAINTSLPFQNIHPITIGKCPKYVKSTKLRLATGLRNVPSIQSRGLFGAIAGFIEGGWTGMVDGWYGYHHQNEQGSGYAADLKSTQNAIDEITNKVNSVIEKMNTQFTAVGKEFNHLEKRIENLNKKVDDGFLDIWTYNAELLVLLENERTLDYHDSNVKNLYEKVRSQLKNNAKEIGNGCFEFYHKCDNTCMESVKNGTYDYPKYSEEAKLNREEIDGVKLESTRIYQILAIYSTVASSLVLVVSLGAISFWMCSNGSLQCRICI', 'H3' : 'QDLPGNDNSTATLCLGHHAVPNGTLVKTITDDQIEVTNATELVQSSSTGKICNNPHRILDGIDCTLIDALLGDPHCDVFQNETWDLFVERSKAFSNCYPYDVPDYASLRSLVASSGTLEFITEGFTWTGVTQNGGSNACKRGPGSGFFSRLNWLTKSGSTYPVLNVTMPNNDNFDKLYIWGIHHPSTNQEQTSLYVQASGRVTVSTRRSQQTIIPNIGSRPWVRGLSSRISIYWTIVKPGDVLVINSNGNLIAPRGYFKMRTGKSSIMRSDAPIDTCISECITPNGSIPNDKPFQNVNKITYGACPKYVKQNTLKLATGMRNVPEKQTRGLFGAIAGFIENGWEGMIDGWYGFRHQNSEGTGQAADLKSTQAAIDQINGKLNRVIEKTNEKFHQIEKEFSEVEGRIQDLEKYVEDTKIDLWSYNAELLVALENQHTIDLTDSEMNKLFEKTRRQLRENAEEMGNGCFKIYHKCDNACIESIRNGTYDHDVYRDEALNNRFQIKGVELKSGYKDWILWISFAISCFLLCVVLLGFIMWACQRGNIRCNICI','B73' : 'DRICTGITSSNSPHVVKTATQGEVNVTGVIPLTTTPTKSHFANLKGTQTRGKLCPNCLNCTDLDVALGRPKCMGTIPSAKASILHEVKPVTSGCFPIMHDRTKIRQLPNLLRGYENIRLSARNVTNAETAPGGPYIVGTSGSCPNVTNGNGFFATMAWAVPKNKTATNPLTVEVPYICTKGEDQITVWGFHSDDETQMVKLYGDSKPQKFTSSANGVTTHYVSQIGGFPNQAEDEGLPQSGRIVVDYMVQKPGKTGTIAYQRGVLLPQKVWCASGRSKVIKGSLPLIGEADCLHEKYGGLNKSKPYYTGEHAKAIGNCPIWVKTPLKLANGTKYRPPAKLLKERGFFGAIAGFLEGGWEGMIAGWHGYTSHGAHGVAVAADLKSTQEAINKITKNLNSLSELEVKNLQRLSGAMDELHNEILELDEKVDDLRADTISSQIELAVLLSNEGIINSEDEHLLALERKLKKMLGPSAVDIGNGCFETKHKCNQTCLDRIAAGTFNAGEFSLPTFDSLNITAASLNDDGLDNHTILLYYSTAASSLAVTLMIAIFIVYMVSRDNVSCSICL', 'PB2' : 'MERIKELRDLMSQSRTREILTKTTVDHMAIIKKYTSGRQEKNPALRMKWMMAMKYPITADKRIIEMIPERNEQGQTLWSKTNDAGSDRVMVSPLAVTWWNRNGPTTSAVHYPKVYKTYFEKVERLKHGTFGPVHFRNQVKIXRRVDINPGHADLSAKEAQDVIMEVVFPNEVGARILTSESQLTITKEKKEELQDCKIAPLMVAYMLERELVRKTRFLPVAGGTSSVYIEVLHLTQGTCWEQMYTPGGEVRNDDVDQSLIIAARNIVRRATVSADPLASLLEMCHSTQIGGTRMVDILRQNPTEEQAVDICKAAMGLRISSSLSFGGFTFKRTSGSSVTKEEEVLTGNLQTLKIRVHEGYEEFTMVGRRATAILRKATRRLIQLIVSGRDEQSIAEAIIVAMVFSQEDCMIKAVRGDLNFVNRANQRLNPMHQLLRHFQKDAKVLFQNWGIEPIDNVMGMIGILPDMTPSTEMSLRGVRVSKMGVDEYSSTERVVVSIDRFLRVRDQRGNVLLSPEEVSETQGIERLTITYSSSMMWEINGPESVLVNTYQWIIRNWETVKIQWSQDPTMLYNKMEFEPFQSLVPKAARGQYSGFVRTLFQQMRDVLGTFDTVQIIKLLPFAAAPPEQSRMQFSSLTVNVRGSGMRILVRGNSPVFNYNKATKRLTVLGKDAGALTEDPDEGTAGVESAVLRGFLILGKEDKRYGPALSINELSNLAKGEKANVLIGQGDVVLVMKRKRDSSILTDSQTATKRIRMAIN', 'PB1' : 'MDVNPTLLFLKVPVQNAISTTFPYTGDPPYSHGTGTGYTMDTVNRTHQYSEKGKWTTNTETGAPQLNPIDGPLPEDNEPSGYAQTDCVLEAMAFLEESHPGIFENSCLETMEIVQQTRVDKLTQGRQTYDWTLNRNQPAATALANTIEIFRSNGLTANESGRLIDFLKDVMESMDKEEMEITTHFQRKRRVRDNMTKKMVTQRTIGKKKQRLNKRSYLIRALTLNTMTKDAERGKLKRRAIATPGMQIRGFVYFVETLARSICEKLEQSGLPVGGNEKKAKLANVVRKMMTNSQDTELSFTITGDNTKWNENQNPRMFLAMITYITRNQPEWFRNVLSIAPIMFSNKMARLGKGYMFESKSMKLRTQIPAEMLASIDLKYFNELTKKKIEKIRPLLIDGAASLSPGMMMGMFNMLSTVLGVSILNLGQKRYTKTTYWWDGLQSSDDFALIVNAPNHEGIQAGVDRFYRTCKLVGINMSKKKSYINRTGTFEFTSFFYRYGFVANFSMELPSFGVSGINESADMSIGVTVIKNNMINNDLGPATAQMALQLFIKDYRYTYRCHRGDTQIQTRRSFELKKLWEQTRSKAGLLVSDGGPNLYNIRNLHIPEVCLKWELMDEDYQGRLCNPLNPFVSHKEIESVNNAVVMPAHGPAKSMEYDAVATTHSWIPKRNRSILNTSQRGILEDEQMYQKCCNLFEKFFPSSSYRRPVGISSMVEAMVSRARIDARIDFESGRIKKEEFAEIMKICSTIEELRRQK', 'PA' : 'MEDFVRQCFNPMIVELAEKAMKEYGEDPKIETNKFAAICTHLEVCFMYSDFHFIDERSESIIVESGDPNALLKHRFEIIEGRDRTMAWTVVNSLCNTTGVEKPKFLPDLYDYKENRFIEIGVTRREVHTYYLEKANKIKSEKTHIHIFSFTGEEMATKADYTIDEESRARIKTRLFTIRQEMASRGLWDSFRQSERGEETIEEKFEITGTMRRLADQSLPPNFSSLENFRAYVDGFEPNGCIEGKLSQMSKEVNARIEPFLKTTPRPLRLPDGPPCSQRSKFLLMDALKLSIEDPSHEGEGIPLYDAIKCMKTFFGWKEPNIVKPHEKGINPNYLLAWKQVLAELQDIENEEKIPKTKNMKKTSQLKWALGENMAPEKVDFEDCKDVSDLTQYNSDEPESRSLASWIQSEFNKACELTDSSWIELDEIGEDVAPIEHIASMRRNYFTAEVSHCRATEYIMKGVYINTALLNASCAAMDDFQLIPMISKCRTKEGRRKTNLYGFIIKGRSHLRNDTDVVNYVSMEFSLTDPRLEPHKWEKYCVLEIGDMLLRTAVGQVSRPMFLYVRTNGTSKIKMKWGMEMRRCLLQSLQQIESMIEAESSVKEKDMTKEFFENKSETWPIGESPKGVEEGSIGKVCRTLLAKSVFNSLYASSQLEGFSAESRKLLLIAQALRDNLEPGTFDLGGLYEAIEECLINDPWVLLNASWFNSFLAHALK'}

    import argparse
    params = argparse.ArgumentParser(description='FoldX python shell script (version 2).')
    params.add_argument('-p', '--pdb', required=True, type=str, help='PDB structure file.')
    params.add_argument('-m', '--mutfile', required=True, type=str, help='FASTA alignment (1st sequence = reference strain and wild-type; Headers will be used as substitution IDs) OR text file with substitutions to perform (every new line in a new substitution; multiple substitutions separated by comma; Can include substitution ID - "ID\\tsubstitution(s)).')
    params.add_argument('-r', '--reference', type=str, default='H3', help='Reference sequence for numbering (default by preset: %(default)s). Presets available: {}. Other acceptable inputs include the sequence itself or FASTA file (will use first sequence entry).'.format(', '.join(preset_reference_sequences.keys())))
    params.add_argument('--repair_pdb', action='store_true', help='Repair PDB input file before performing stability analysis.')
    params.add_argument('--chains', default='1', choices=map(str, [1,2,3]), help='Chain options (default %(default)s). 1: Single protein broken into multiple chains, 2: Multiple chains of the SAME protein, 3: Multiple chains of DIFFERENT proteins.')
    params.add_argument('--sequential', action='store_true', help='Change analysis mode to sequential mode.')
    params.add_argument('--save_mt_pdb', action='store_true', help='Save mutant pdb files.')
    params.add_argument('--outfname', type=str, help='Optional output filename.')
    # foldx parameters
    params.add_argument('--runs', type=int, help='No. of replicates (Default: %(default)s)', default=5)
    params.add_argument('--pH', type=float, help='pH (Default: %(default)s)', default=7.0)
    params.add_argument('--temp', type=int, help='Temperature in Kelvins (Default: %(default)sK)', default=298)
    params.add_argument('--ion', type=float, help='Ionic strength in M (Default: %(default)sM)', default=0.05)
    params.add_argument('--move_neighbors', action='store_true', help='Move neighboring atoms.')
    params = params.parse_args()

    # parse pdb
    import re
    writepdbBIN, pdbfname = re.search('(/*)([^\/]+)\.pdb$', params.pdb).group(1, 2)

    # if pdb file is not in current working directory, copy it over
    from os.path import expanduser
    import subprocess
    if writepdbBIN != '':
        cmd = ['cp', expanduser(params.pdb), '{}.pdb'.format(pdbfname)]
        subprocess.call(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

    # std out file for all subprocesses
    stdout_file = open('foldx-output_verbose_{}.txt'.format(pdbfname), 'w')

    # Repair pdb if required
    import sys, platform
    # get home directory
    home = expanduser('~')
    # check platform
    system = 'mac' if platform.system() == 'Darwin' else 'linux'
    if params.repair_pdb:
        print ('\nRepairing pdb...')
        cmd = [expanduser('~/Dropbox/softwares/foldx_{}/foldx'.format(system)), '--command=RepairPDB',
               '--pdb={}.pdb'.format(pdbfname),
               '--rotabaseLocation={}/Dropbox/softwares/foldx_{}/rotabase.txt'.format(home, system)]
        subprocess.call(cmd, stderr=subprocess.PIPE, stdout=stdout_file)

        pdbfname = '{}_Repair'.format(pdbfname)
        print ('...completed.')

    # parse pdb file
    try:
        chain_to_seqindex_to_pos_res, chain_to_sequence = parsepdb('{}.pdb'.format(pdbfname))
    except:
        sys.exit('\nERROR: Invalid PDB file.\n')

    # output filename
    if not params.outfname:
        params.outfname = 'foldx-output{}_{}_{}.txt'.format('-sequential' if params.sequential else '', pdbfname, re.sub('([^/]+/|\.fa$|\.txt$)', '', params.mutfile))
    print ('\noutput file...{}'.format(params.outfname))

    # parse mutfile
    import itertools
    substitution_list = []
    try:
        substitution_list_fdat, substitution_id = parsefasta(params.mutfile, 1)
        print ('\ninput mutfile..FASTA.')
        # first sequence is WT and reference for numbering
        params.reference = substitution_list_fdat[substitution_id[0]]

        id_to_remove = []
        for _ in xrange(1, len(substitution_id), 1):
            curr_sequence = substitution_list_fdat[substitution_id[_]]
            if params.sequential:
                prev_sequence = substitution_list_fdat[substitution_id[_-1]]
            else:
                prev_sequence = params.reference

            substitutions = tuple(['{}{}{}'.format(prev_sequence[i], i+1, curr_sequence[i]) for i in xrange(len(params.reference)) if prev_sequence[i] != curr_sequence[i]])
            if len(substitutions) > 0:
                substitution_list.append(substitutions)
            else:
                # no substitution - remove id later
                id_to_remove.append(_)

        # remove id with no substitutions
        if len(id_to_remove) > 0:
            for id in id_to_remove[::-1]:
                del substitution_id[id]

        # remove first reference id
        substitution_id.pop(0)

    except:
        substitution_id = []
        print ('\ninput mutfile..TEXT.')
        try:
            fhandle = filter(None, [_.strip() for _ in open(params.mutfile, 'rU')])
            for i, line in enumerate(fhandle):
                try:
                    id, line = line.split('\t')
                except:
                    id = 'SUBID_{:03d}'.format(i+1)

                try:
                    substitutions = tuple(re.findall('[A-Z]\d+[A-Z]', line))
                except:
                    sys.exit()

                substitution_list.append(substitutions)
                substitution_id.append(id)
        except:
            sys.exit('\nERROR: Invalid input mutfile.\n')

    # remove illegal characters from ID
    substitution_id = [re.sub('\.', 'p', _) for _ in substitution_id]

    # parse for reference sequence for numbering
    try:
        reference_sequence = preset_reference_sequences[params.reference]
        print ('preset reference sequence...{}'.format(params.reference))
    except:
        try:
            reference_sequence = re.search('^[rhkdestnqcgpavilmfyw]+$', params.reference, re.I).group().upper()
            print ('custom reference sequence given...')
        except:
            try:
                ref_fdat, ref_headers = parsefasta(params.reference, 1)
                reference_sequence = re.search('^[rhkdestnqcgpavilmfyw]+$', ref_fdat[ref_headers[0]], re.I).group().upper()
            except:
                sys.exit('\nERROR: Unable to parse reference seqeunce.\n')

    # match chain to reference proteins
    print ('...aligning structure and reference sequences...')
    refpos_to_chain_structurepos_res = {}
    # if the same molecule is broken into chains - align all chains to the same protein
    if params.chains in map(str, [1,2]):
        for chain, sequence in chain_to_sequence.items():
            with open('temp_{}.fa'.format(chain), 'w') as output:
                output.write('>reference\n{}\n>{}\n{}\n'.format(reference_sequence, chain, sequence))
            subprocess.call('mafft --maxiterate 1000 --globalpair temp_{chain}.fa > temp_{chain}.mafft.fa'.format(chain=chain), shell=True, stderr=subprocess.STDOUT, stdout=stdout_file)

            for refpos, (structure_pos, structure_res) in get_ref_to_structural_pos(chain, parsefasta('temp_{}.mafft.fa'.format(chain))).items():
                if params.chains == '1' and refpos in refpos_to_chain_structurepos_res:
                    sys.exit('\nERROR: Position {} found in > 1 chain. Check alignment.\n'.format(refpos))

                try:
                    refpos_to_chain_structurepos_res[refpos][chain] = (structure_pos, structure_res)
                except:
                    refpos_to_chain_structurepos_res[refpos] = {chain:(structure_pos, structure_res)}

    # if every chain is a different protein - use chains that are best-scoring
    else:
        chain_to_gapcounts_and_fdat = {}
        for chain, sequence in chain_to_sequence.items():
            with open('temp_{}.fa'.format(chain), 'w') as output:
                output.write('>reference\n{}\n>{}\n{}\n'.format(reference_sequence, chain, sequence))
            subprocess.call('mafft --maxiterate 1000 --globalpair temp_{chain}.fa > temp_{chain}.mafft.fa'.format(chain=chain), shell=True, stderr=subprocess.STDOUT, stdout=stdout_file)

            chain_fdat = parsefasta('temp_{}.mafft.fa'.format(chain))
            chain_to_gapcounts_and_fdat[chain] = (chain_fdat[chain].count('-'), chain_fdat)

        # assign best (lowest gaps) chains to reference polymerase proteins
        best_score = min([_[0] for _ in chain_to_gapcounts_and_fdat.values()])

        for chain, (score, chain_fdat) in chain_to_gapcounts_and_fdat.items():
            if score == best_score:
                for refpos, (structure_pos, structure_res) in get_ref_to_structural_pos(chain, chain_fdat).items():
                    try:
                        refpos_to_chain_structurepos_res[refpos][chain] = (structure_pos, structure_res)
                    except:
                        refpos_to_chain_structurepos_res[refpos] = {chain:(structure_pos, structure_res)}

    # write output headers
    with open(params.outfname, 'w') as output:
        output.write('{}SUB_ID\tMain(*)\tWT<>Struc-Res(*)\tSubstitution(s)\tWT-mean\tWT-SD\tMT-mean\tMT-SD\tddG-mean\tddG-SD\n'.format('Step\t' if params.sequential else ''))

    # perform FoldX stability analyses
    print ('\n...performing FoldX analyses...')
    import numpy as np

    if params.sequential:
        print ('...sequential anlaysis...')
        ### --- sequential --- ###
        refpos_to_residue_order = {}
        id_to_no_foldx = {}
        WT_substitutions_to_make = {}
        individual_list = {}

        for id_index, id in enumerate(substitution_id):
            curr_substitutions = substitution_list[id_index]

            for sub in curr_substitutions:
                WT, refpos, MT = re.search('([A-Z])(\d+)([A-Z])', sub).group(1,2,3)
                refpos = int(refpos)

                # check for WT substitution to make
                if refpos in refpos_to_residue_order:
                    # double check that substitutions are carried out in sequential order
                    if WT == refpos_to_residue_order[refpos][-1]:
                        refpos_to_residue_order[refpos].append(MT)
                    else:
                        sys.exit('\nERROR: WT residue ({}) in position {} for {} is not in line with previous MT ({}, {})\n'.format(WT, refpos, id, refpos_to_residue_order[refpos][-1], substitution_id[id_index-1]))

                else:
                    refpos_to_residue_order[refpos] = [WT, MT]
                    if refpos in refpos_to_chain_structurepos_res:
                        for chain, (structure_pos, structure_res) in refpos_to_chain_structurepos_res[refpos].items():
                            # check for WT substitutions to make
                            if structure_res != WT:
                                WT_sub = '{}{}{}'.format(structure_res, refpos, WT)
                                try:
                                    WT_substitutions_to_make[WT_sub].append('{}{}{}{}'.format(structure_res, chain, structure_pos, WT))
                                except:
                                    WT_substitutions_to_make[WT_sub] = ['{}{}{}{}'.format(structure_res, chain, structure_pos, WT)]

                # parse for individual list
                if refpos in refpos_to_chain_structurepos_res:
                    for chain, (structure_pos, structure_res) in refpos_to_chain_structurepos_res[refpos].items():
                        # append substitution to individual list
                        try:
                            individual_list[id][sub].append('{}{}{}{}'.format(WT, chain, structure_pos, MT))
                        except:
                            try:
                                individual_list[id][sub] = ['{}{}{}{}'.format(WT, chain, structure_pos, MT)]
                            except:
                                individual_list[id] = {sub:['{}{}{}{}'.format(WT, chain, structure_pos, MT)]}
                else:
                    # position not found in structure
                    if len(curr_substitutions) == 1:
                        # if there is only 1 substitution, print to output file and set perform_foldx to zero
                        id_to_no_foldx[id] = 0
                    else:
                        try:
                            id_to_no_foldx[id].append(sub)
                        except:
                            id_to_no_foldx[id] = [sub]

        # id_index to start after first substitution (can be WT sub or first sub)
        start_id_index = -1

        # first make WT substitutions if any
        if len(WT_substitutions_to_make) > 0:
            result = foldx([x for y in WT_substitutions_to_make.values() for x in y], params.runs, pdbfname, params.ion, params.pH, params.temp, 1)
            WT_values = [result[r]['WT'] for r in xrange(1, params.runs+1, 1)]
            MT_values = [result[r]['MT'] for r in xrange(1, params.runs+1, 1)]
            ddG_values = [result[r]['MT']-result[r]['WT'] for r in xrange(1, params.runs+1, 1)]

            # write results to output
            with open(params.outfname, 'a') as output:
                output.write('-1\tWT-SUB\t\t*\t{}\t{}\n'.format(','.join(WT_substitutions_to_make.keys()), '\t'.join(map(str, [np.mean(WT_values), np.std(WT_values, ddof=1), np.mean(MT_values), np.std(MT_values, ddof=1), np.mean(ddG_values), np.std(ddG_values, ddof=1)]))))

            # save mt structures
            for r in xrange(params.runs):
                cmd = 'mv {}_1_{}.pdb {}_WT-SUB_R{}.pdb'.format(pdbfname, r, pdbfname, r+1)
                subprocess.call(cmd, shell=True, stderr=subprocess.STDOUT, stdout=stdout_file)
            prev_structure_fname = '{}_WT-SUB'.format(pdbfname)

        # if there are no WT substitutions, perform foldx for first index-ed substitution(s)
        else:
            for id_index, id in enumerate(substitution_id):
                curr_substitutions = substitution_list[id_index]

                if id in id_to_no_foldx:
                    if id_to_no_foldx[id] == 0:
                        # continue if no substitution is being made in current step
                        with open(params.outfname, 'a') as output:
                            output.write('{}\t{}\t*\t\t{}\tPosition {} not found in given structure.\n'.format(id_index, id, curr_substitutions[0], re.search('\d+', curr_substitutions[0]).group()))
                        continue
                    else:
                        # if individual substitution's position cannot be found
                        for indv_id in id_to_no_foldx[id]:
                            with open(params.outfname, 'a') as output:
                                output.write('{}\t{}\t\t\t{}\tPosition {} not found in given structure.\n'.format(id_index, id, indv_id, re.search('\d+', indv_id).group()))

                for indv_id, indv_sub in individual_list[id].items():
                    # perform each individual substitution first
                    result = foldx(indv_sub, params.runs, pdbfname, params.ion, params.pH, params.temp, 1)
                    WT_values = [result[r]['WT'] for r in xrange(1, params.runs+1, 1)]
                    MT_values = [result[r]['MT'] for r in xrange(1, params.runs+1, 1)]
                    ddG_values = [result[r]['MT']-result[r]['WT'] for r in xrange(1, params.runs+1, 1)]

                    # write results to output
                    with open(params.outfname, 'a') as output:
                        # if there is only 1 substitution
                        if len(curr_substitutions) == 1:
                            output.write('{}\t{}\t*\t\t{}\t{}\n'.format(id_index, id, indv_id, '\t'.join(map(str, [np.mean(WT_values), np.std(WT_values, ddof=1), np.mean(MT_values), np.std(MT_values, ddof=1), np.mean(ddG_values), np.std(ddG_values, ddof=1)]))))
                        else:
                            output.write('{}\t{}\t\t\t{}\t{}\n'.format(id_index, id, indv_id, '\t'.join(map(str, [np.mean(WT_values), np.std(WT_values, ddof=1), np.mean(MT_values), np.std(MT_values, ddof=1), np.mean(ddG_values), np.std(ddG_values, ddof=1)]))))

                            # also remove mt pdb files if we are not working on the main substitution already
                            cmd = 'rm {}_1_*.pdb'.format(pdbfname)
                            subprocess.call(cmd, shell=True, stderr=subprocess.STDOUT, stdout=stdout_file)

                # if there are more than 1 substitution component
                if len(curr_substitutions) > 1:
                    # perform for all substitutions
                    result = foldx([x for y in individual_list[id].values() for x in y], params.runs, pdbfname, params.ion, params.pH, params.temp, 1)
                    WT_values = [result[r]['WT'] for r in xrange(1, params.runs+1, 1)]
                    MT_values = [result[r]['MT'] for r in xrange(1, params.runs+1, 1)]
                    ddG_values = [result[r]['MT']-result[r]['WT'] for r in xrange(1, params.runs+1, 1)]

                    # write results to output
                    with open(params.outfname, 'a') as output:
                        output.write('{}\t{}\t*\t\t{}\t{}\n'.format(id_index, id, ','.join(curr_substitutions), '\t'.join(map(str, [np.mean(WT_values), np.std(WT_values, ddof=1), np.mean(MT_values), np.std(MT_values, ddof=1), np.mean(ddG_values), np.std(ddG_values, ddof=1)]))))

                # save mt structures
                for r in xrange(params.runs):
                    cmd = 'mv {}_1_{}.pdb {}_{}_R{}.pdb'.format(pdbfname, r, pdbfname, id, r+1)
                    subprocess.call(cmd, shell=True, stderr=subprocess.STDOUT, stdout=stdout_file)
                prev_structure_fname = '{}_{}'.format(pdbfname, id)

                # update start_id_index
                start_id_index = id_index
                break

        # continue with the rest of the substitutions
        for id_index, id in enumerate(substitution_id):
            # only start after start_id_index
            if id_index <= start_id_index:
                continue

            curr_substitutions = substitution_list[id_index]

            if id in id_to_no_foldx:
                if id_to_no_foldx[id] == 0:
                    # continue if no substitution is being made in current step
                    with open(params.outfname, 'a') as output:
                        output.write('{}\t{}\t*\t\t{}\tPosition {} not found in given structure.\n'.format(id_index, id, curr_substitutions[0], re.search('\d+', curr_substitutions[0]).group()))
                    continue
                else:
                    # if individual substitution's position cannot be found
                    for indv_id in id_to_no_foldx[id]:
                        with open(params.outfname, 'a') as output:
                            output.write('{}\t{}\t\t\t{}\tPosition {} not found in given structure.\n'.format(id_index, id, indv_id, re.search('\d+', indv_id).group()))

            for indv_id, indv_sub in individual_list[id].items():
                # perform each individual substitution first
                result = {}
                for r in xrange(1, params.runs+1, 1):
                    pdbfname_to_use = '{}_R{}'.format(prev_structure_fname, r)
                    result[r] = foldx(indv_sub, 1, pdbfname_to_use, params.ion, params.pH, params.temp, 1)[1]

                    if len(curr_substitutions) == 1:
                        # save as next run WT pdb
                        cmd = 'mv {}_1.pdb {}_{}_R{}.pdb'.format(pdbfname_to_use, pdbfname, id, r)
                        subprocess.call(cmd, shell=True, stderr=subprocess.STDOUT, stdout=stdout_file)
                    else:
                        # also remove mt pdb files if we are not working on the main substitution already
                        cmd = 'rm {}_1.pdb'.format(pdbfname_to_use)
                        subprocess.call(cmd, shell=True, stderr=subprocess.STDOUT, stdout=stdout_file)

                WT_values = [result[r]['WT'] for r in xrange(1, params.runs+1, 1)]
                MT_values = [result[r]['MT'] for r in xrange(1, params.runs+1, 1)]
                ddG_values = [result[r]['MT']-result[r]['WT'] for r in xrange(1, params.runs+1, 1)]

                # write results to output
                with open(params.outfname, 'a') as output:
                    # if there is only 1 substitution component
                    if len(curr_substitutions) == 1:
                        output.write('{}\t{}\t*\t\t{}\t{}\n'.format(id_index, id, indv_id, '\t'.join(map(str, [np.mean(WT_values), np.std(WT_values, ddof=1), np.mean(MT_values), np.std(MT_values, ddof=1), np.mean(ddG_values), np.std(ddG_values, ddof=1)]))))
                    else:
                        output.write('{}\t{}\t\t\t{}\t{}\n'.format(id_index, id, indv_id, '\t'.join(map(str, [np.mean(WT_values), np.std(WT_values, ddof=1), np.mean(MT_values), np.std(MT_values, ddof=1), np.mean(ddG_values), np.std(ddG_values, ddof=1)]))))

            # if there are more than 1 substitution component
            if len(curr_substitutions) > 1:
                result = {}
                for r in xrange(1, params.runs+1, 1):
                    pdbfname_to_use = '{}_R{}'.format(prev_structure_fname, r)
                    result[r] = foldx([x for y in individual_list[id].values() for x in y], 1, pdbfname_to_use, params.ion, params.pH, params.temp, 1)[1]

                    # save as next run WT pdb
                    cmd = 'mv {}_1.pdb {}_{}_R{}.pdb'.format(pdbfname_to_use, pdbfname, id, r)
                    subprocess.call(cmd, shell=True, stderr=subprocess.STDOUT, stdout=stdout_file)

                WT_values = [result[r]['WT'] for r in xrange(1, params.runs+1, 1)]
                MT_values = [result[r]['MT'] for r in xrange(1, params.runs+1, 1)]
                ddG_values = [result[r]['MT']-result[r]['WT'] for r in xrange(1, params.runs+1, 1)]

                # write results to output
                with open(params.outfname, 'a') as output:
                    output.write('{}\t{}\t*\t\t{}\t{}\n'.format(id_index, id, ','.join(curr_substitutions), '\t'.join(map(str, [np.mean(WT_values), np.std(WT_values, ddof=1), np.mean(MT_values), np.std(MT_values, ddof=1), np.mean(ddG_values), np.std(ddG_values, ddof=1)]))))

            # update prev_structure_fname
            prev_structure_fname = '{}_{}'.format(pdbfname, id)

    else:
        ### --- normal --- ###
        print ('...normal analysis...')
        for id_index, id in enumerate(substitution_id):
            curr_substitutions = substitution_list[id_index]

            perform_foldx = 1

            WT_substitutions_to_make = {}
            individual_list = {}

            for sub in curr_substitutions:
                WT, refpos, MT = re.search('([A-Z])(\d+)([A-Z])', sub).group(1,2,3)
                refpos = int(refpos)

                if refpos in refpos_to_chain_structurepos_res:
                    for chain, (structure_pos, structure_res) in refpos_to_chain_structurepos_res[refpos].items():
                        # check for WT substitutions to make
                        if structure_res != WT:
                            WT_sub = '{}{}{}'.format(structure_res, refpos, WT)
                            try:
                                WT_substitutions_to_make[WT_sub].append('{}{}{}{}'.format(structure_res, chain, structure_pos, WT))
                            except:
                                WT_substitutions_to_make[WT_sub] = ['{}{}{}{}'.format(structure_res, chain, structure_pos, WT)]

                        # append substitution to individual list
                        if len(curr_substitutions) == 1:
                            # if there is only 1 substitution
                            try:
                                individual_list[id].append('{}{}{}{}'.format(WT, chain, structure_pos, MT))
                            except:
                                individual_list[id] = ['{}{}{}{}'.format(WT, chain, structure_pos, MT)]
                        else:
                            try:
                                individual_list[sub].append('{}{}{}{}'.format(WT, chain, structure_pos, MT))
                            except:
                                individual_list[sub] = ['{}{}{}{}'.format(WT, chain, structure_pos, MT)]
                else:
                    # position not found in structure
                    if len(curr_substitutions) == 1:
                        # if there is only 1 substitution, print to output file and set perform_foldx to zero
                        with open(params.outfname, 'a') as output:
                            output.write('{}\t*\t\t{}\tPosition {} not found in given structure.\n'.format(id, curr_substitutions[0], refpos))
                        perform_foldx = 0

                    else:
                        # print to output file
                        with open(params.outfname, 'a') as output:
                            output.write('{}\t\t\t{}\tPosition {} not found in given structure.\n'.format(id, sub, refpos))

            # no need to perform foldx if there is only 1 substitution and position is not found in structure
            if perform_foldx == 0:
                continue

            # consolidate all substitutions
            if len(curr_substitutions) > 1:
                individual_list[id] = [x for y in individual_list.values() for x in y]

            # make WT substitutions first if any
            if len(WT_substitutions_to_make) > 0:
                result = foldx([x for y in WT_substitutions_to_make.values() for x in y], params.runs, pdbfname, params.ion, params.pH, params.temp, 1)
                WT_values = [result[r]['WT'] for r in xrange(1, params.runs+1, 1)]
                MT_values = [result[r]['MT'] for r in xrange(1, params.runs+1, 1)]
                ddG_values = [result[r]['MT']-result[r]['WT'] for r in xrange(1, params.runs+1, 1)]

                # write results to output
                with open(params.outfname, 'a') as output:
                    output.write('{}\t\t*\t{}\t{}\n'.format(id, ','.join(WT_substitutions_to_make.keys()), '\t'.join(map(str, [np.mean(WT_values), np.std(WT_values, ddof=1), np.mean(MT_values), np.std(MT_values, ddof=1), np.mean(ddG_values), np.std(ddG_values, ddof=1)]))))

            for indv_id, indv_sub in individual_list.items():
                # if we have made any WT substitution, use that structure
                if len(WT_substitutions_to_make) > 0:
                    result = {}
                    # analyse each substitution by single runs
                    for r in xrange(params.runs):
                        pdbfname_to_use = '{}_1_{}'.format(pdbfname, r)
                        result[r+1] = foldx(indv_sub, 1, pdbfname_to_use, params.ion, params.pH, params.temp, 1)[1]

                        if params.save_mt_pdb and indv_id == id:
                            # save final pdb
                            cmd = 'mv {}_1.pdb {}_{}_R{}.pdb'.format(pdbfname_to_use, pdbfname, id, r+1)
                            subprocess.call(cmd, shell=True, stderr=subprocess.STDOUT, stdout=stdout_file)
                        else:
                            cmd = 'rm {}_1.pdb'.format(pdbfname_to_use)
                            subprocess.call(cmd, shell=True, stderr=subprocess.STDOUT, stdout=stdout_file)
                else:
                    pdbfname_to_use = pdbfname
                    result = foldx(indv_sub, params.runs, pdbfname_to_use, params.ion, params.pH, params.temp, 1)
                    
                    if params.save_mt_pdb and indv_id == id:
                        # save final pdb
                        for r in xrange(params.runs):
                            cmd = 'mv {}_1_{}.pdb {}_{}_R{}.pdb'.format(pdbfname_to_use, r, pdbfname, id, r+1)
                            subprocess.call(cmd, shell=True, stderr=subprocess.STDOUT, stdout=stdout_file)
                    else:
                        cmd = 'rm {}_1_*.pdb'.format(pdbfname_to_use)
                        subprocess.call(cmd, shell=True, stderr=subprocess.STDOUT, stdout=stdout_file)

                # consolidate results
                WT_values = [result[r]['WT'] for r in xrange(1, params.runs+1, 1)]
                MT_values = [result[r]['MT'] for r in xrange(1, params.runs+1, 1)]
                ddG_values = [result[r]['MT']-result[r]['WT'] for r in xrange(1, params.runs+1, 1)]

                # write results to output
                with open(params.outfname, 'a') as output:
                    output.write('{}\t{}\t\t{}\t{}\n'.format(id, '*' if indv_id == id else '', ','.join(curr_substitutions) if indv_id == id else indv_id, '\t'.join(map(str, [np.mean(WT_values), np.std(WT_values, ddof=1), np.mean(MT_values), np.std(MT_values, ddof=1), np.mean(ddG_values), np.std(ddG_values, ddof=1)]))))

            if len(WT_substitutions_to_make) > 0:
                # remove the mutated-to-WT pdb
                cmd = 'rm {}_1_*.pdb'.format(pdbfname)
                subprocess.call(cmd, shell=True, stderr=subprocess.STDOUT, stdout=stdout_file)

    # remove all mt files if called to be saved
    if params.sequential and not params.save_mt_pdb:
        cmd = 'rm {}_*_R*.pdb'.format(pdbfname)
        subprocess.call(cmd, shell=True, stderr=subprocess.STDOUT, stdout=stdout_file)

    # close stdout output file
    stdout_file.close()

    print ('\n...done.\n')
    exit(0)