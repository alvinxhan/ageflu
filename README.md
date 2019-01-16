# Custom scripts and files for Han et al. (2018)

## Citation

Han A. X., Maurer-Stroh S., Russell C. A.. (2018). Individual immune selection pressure has limited impact on seasonal influenza virus evolution. _Nature Ecology and Evolution_.

## Jupyter notebooks 

* [**age_analysis_notebook.ipyb**](https://github.com/alvinxhan/ageflu/blob/master/age_analysis_notebook.ipynb): Demonstrates how evolutionarily closest sequence pairs were parsed from a phylogenetic tree and perform age-associated analyses on the amino acid substitutions inferred from these closest pairs. 

* [**pair_distribution_analyses.ipyb**](https://github.com/alvinxhan/ageflu/blob/master/pair_distribution_analyses.ipynb): Demonstrates how pair distribution analyses were performed. 

## Files 
You will find the data and result files for all seasonal influenza subtypes/lineages analysed, each with their own subfolder containing:  
_Data files_
* {subtype}-HA-nuc_rooted.nwk - Rooted phylogenetic tree inferred using our nested phylogenetic inference approach. 
* {subtype}-HA-nuc_cd-hit-non-redundant.clstr - Output file of CD-HIT clustering viruses with identical HA nucleotide sequences. 
* {subtype}-HA-nuc.fa - Curated FASTA HA nucleotide sequences. - **We are unable to provide the sequence alignment as GISAID disallows releasing genetic sequences to other external public domains. However, we have provided the isolate identifiers/accession numbers as well as headers to the sequences we have used in Sequence_headers_ALLSUBTYPES.csv and Acknowledge_Table_GISAID.csv.**  

_Result files_

From [ageflu_getpairs.py](https://github.com/alvinxhan/ageflu/blob/master/scripts/ageflu_getpairs.py) in [scripts](https://github.com/alvinxhan/ageflu/blob/master/scripts/) OR if you go through [**age_analysis_notebook.ipyb**](https://github.com/alvinxhan/ageflu/blob/master/age_analysis_notebook.ipynb):  
  * ageflu_evol-closest-pairs.{child-min-age}C{child-max-age}\_{adult-min-age}A{adult-max-age}\_PD0.007_MM5_CL3_{subtype}-HA-nuc_rooted.nwk.txt - Tab-delimited file containing all closely-related virus pairs based on the age limits stipulated in the file title. 
  
From [ageflu_analysis.py](https://github.com/alvinxhan/ageflu/blob/master/scripts/ageflu_analysis.py) in [scripts](https://github.com/alvinxhan/ageflu/blob/master/scripts/) OR if you go through [**age_analysis_notebook.ipyb**](https://github.com/alvinxhan/ageflu/blob/master/age_analysis_notebook.ipynb):
  * ageflu_pairs-distribution.{child-min-age}C{child-max-age}\_{adult-min-age}A{adult-max-age}\_PD0.007_MM5_CL3_{subtype}-HA-nuc_rooted.nwk_SP{same-passage-binary}.txt - Distribution of pairs and corresponding association analyses to determine if there were significant differences between end-in-adult and end-in-child pair types. 
  * ageflu_age-associated-sites.{child-min-age}C{child-max-age}\_{adult-min-age}A{adult-max-age}\_PD0.007_MM5_CL3_{subtype}-HA-nuc_rooted.nwk_SP{same-passage-binary}.txt - Association analyses for every HA amino acid position to determine if there were particular sites that were significantly associated with age. 
  
## Scripts 
You will find custom Python 2 and R scripts used in our analyses here: 
* ageflu_getpairs.py - python script to parse for evolutionarily closest pairs. 
* ageflu_analysis.py - python script to perform association analyses on evolutionarily closest pairs parsed by ageflu_getpairs.py. 
* physicochemical_properties_association_v2.py and foldx_stability_master_v2.py - these two python scripts (the former being a wrapper) were used alongside FoldX to determine if amino acid substitutions with physicochemical property changes were significantly more probable to be resulted from end-in-adult or end-in-child pairs. If you are interested in how to use these scripts, do drop me an email at hanxc@bii.a-star.edu.sg. 
  

* bayesian_binomial_MCMC.R - A simple Bayesian MCMC used to analyse the CD-HIT clstr files to determine if sequences identical at the nucleotide level were equally likely to be found in both adults and children.
* ageflu_gtest.R, ageflu_barnard.R, agresti_and_min.R - These R scripts are the actual hypothesis tests/statistical analysis used in our association analyses. They are embedded within ageflu_analysis.py.

