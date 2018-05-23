# Custom scripts and files for Han et al. (2018)

Citation

Alvin Xiaochuan Han, Sebastian Maurer-Stroh, Colin A. Russell. Individual immune selection pressure has limited impact on seasonal influenza virus evolution.

## Usage instructions
1. Parse for evolutionarily closely-related virus pairs from the inferred phylogenetic tree (rooted, in NEWICK format) using ageflu_getpairs.py. Additional inputs include: (i) FASTA alignment of all sequences (including identical nucleotide sequences), (ii) CD-HIT cluster file of non-redundant sequence clusters. 
  * The age limits of children and adults can be changed using the '--child' and '--adult' arguments. 
  * Note that the FASTA alignment containing reference sequences for HA-numbering ('H1pdm09_H3_FluB_NumberingRef.fa') must be in the same folder as the ageflu_getpairs.py script. 
  * Tab-delimited output: ageflu_evol-closest-pairs.\*.txt

2. Using output from (1) and FASTA alignment of all sequences, use ageflu_analysis.py to get pairs distributions and perform all association analyses as detailed in the paper. The script will generate the following outputs: 
  * ageflu_pairs-distribution.\*.txt:  Pair distributions and substitution frequency analyses.
  * ageflu_age-associated-sites.\*.txt: Sites with substitution frequencies associated with age and their corresponding association analyses.
