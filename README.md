[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.10960421.svg)](https://doi.org/10.5281/zenodo.10960421)

# Supporting data and code for: *Myzus persicae* resistance to neonicotinoids - unravelling the contribution of different mechanisms to phenotype.
*This repository contains the R code used for the data analyses and production of the figures of the related article*

![alt text](https://db3pap005files.storage.live.com/y4mDlIdRciVaiIk5fs2Dl1bMvduiUuhuSt5MaPMCfoKs6iG_Wv1pzwiD5bzqd8Cc0khCc3R2vKX_Dy0G7SA_ZAKBkK7Z6ErM82rmN0PBrGVhi5J_X1MgkOaWnPbpG-9XsdJlQBP_EH7NkxoujphX2IPsmmx_-znhZq-RUVToP9UaiIKUqIxeI-bNseRTtwMZXvy?width=1584&height=588&cropmode=none)


## Context
Insecticide resistance can lead to the repeated application of treatments that are both ineffective in pest control and harmful for the environment. In order to establish sound strategies to slow down or limit these resistance phenomena, precise knowledge of the underlying mechanisms, their relative contribution and their potential interaction is essential. *Myzus persicae* is a major pest aphid, capable of infesting a wide range of crops, with significant economic impact. *M. persicae* is known to have evolved resistances to insecticides of different families, including the abundantly used neonicotinoids. *M. persicae* resistance to neonicotinoids has previously been described as being due to two main mechanisms: a P450 overexpression metabolic resistance and a target-site mutation, R81T, but their respective contribution to resistant phenotypes remains unclear. In this study we combined extensive bioassays with synergist on numerous clones, with gene copy number and expression quantification of two key P450 enzymes (*CYP6CY3* and *CYP6CY4*) to explain the observed phenotypes and assess the relative contribution of metabolic and target-site mechanisms to neonicotinoid resistance in *M. persicae*.


## Datasets
There are 2 datasets used in this study. The files can be found in the "data" folder. 

+ **bioassayRawData.txt:** the first dataset contains the raw data of the results of the bioassays. 
  + *ana_id*: analysis ID
  + *ech_id*: sample ID
  + *prog_id*: program ID
  + *bioagr_id*: pest species name
  + *hote_id*: host name (in French)
  + *pest_sa_id*: name of the pesticide active substance
  + *pest_pc_id*: name of the plant protection product
  + *synerg_id*: name of the synergist used (=AUCUN if none have been used)
  + *tech_id*: name of the bioassay method used (in French)
  + *tps_expo*: duration of exposure
  + *bioagr_class_id*: trait or category evaluated
  + *bioagr_class_val*: value of the trait or category
  + *sex*: for species for which the sex as an effect on resistance
  + *dat_test*: date of the bioassay
  + *rep_test*: replicate number (if any)
  + *lect_id*: reading ID of the bioassay (if there is multiple readings)
  + *test_echec*: was the bioassay successful or not (0=failed, 1=succes)
  + *dose*: active substance dose (mg/L)
  + *nb_vi*: number of individual alive
  + *nb_mb*: number of individual moribund
  + *nb_mt*: number of dead individual
  + *nb_mtot*: total number of individual dead or moribund (=nb_mb+nb_mt)

+ **summaData.txt:** the second dataset contains the phenotypic and genotypic information for the tested clones. Each line correspond to a clone. 
  + *clone-ID*: identification number of the clone
  + *genetic-group*: the Neutral genetic cluster information. Either primary host ('PeaHost') or secondary host ('SecHost')
  + *nAChR-81*: genotype at the R81T loci. 'RR' for homozygous sensitive, 'TT' for homozygous resistant and 'RT' for heterozygous
  + *LC50*: thiaclopride median lethal dose for the clone
  + *LC50-PBO*: thiaclopride + PBO median lethal dose for the clone
  + *CY3_CN*: CYP6CY3 gene copy number evaluated with quantitative PCR
  + *CY4_CN*: CYP6CY4 gene copy number evaluated with quantitative PCR
  + *CY23_CN*: CYP6CY23 gene copy number evaluated with quantitative PCR
  + *CY3_EXP*: CYP6CY3 expression level evaluated with RT-quantitative PCR
  + *CY4_EXP*: CYP6CY4 expression level evaluated with RT-quantitative PCR
  + *CY23_EXP*: CYP6CY23 expression level evaluated with RT-quantitative PCR
  + *CY3_SE_CN*: Standard error of CYP6CY3 gene copy number evaluated with RT-quantitative PCR
  + *CY4_SE_CN*: Standard error of CYP6CY4 gene copy number evaluated with RT-quantitative PCR
  + *CY23_SE_CN*: Standard error of CYP6CY23 gene copy number evaluated with RT-quantitative PCR
  + *CY3_SE_EXP*: Standard error of CYP6CY3 expression level evaluated with RT-quantitative PCR
  + *CY4_SE_EXP*: Standard error of CYP6CY4 expression level evaluated with RT-quantitative PCR
  + *CY23_SE_EXP*: Standard error of CYP6CY23 expression level evaluated with RT-quantitative PCR


## R scripts
+ **MyMetab_load.R:** the script to load the different datasets and packages in the environment. This must be run before any other scripts. 
+ **MyMetab_drc.R:** this script analyses the raw data of the bioassay experiment (see 'bioassayRawData.txt' data file) to compute the thiaclorpide LD50 for each clone with or without the addition of PBO. It produce two files: a pdf files with a page for each clone summarizing the results of the dose-response analyses ('figure_by_clone.pdf') and a table which gives for each clone: the LD50, the standard error and the mean number of individual tested per dose with or whithout PBO, as well as the results of test comparing the LD50 values with and whithout PBO ('results_bioassay.txt').
+ **MyMetab_glm.R:** this script described the building of the glm to model the LD50 as a function of variable of interests. The first step consists in investigating the correlation between P450 quantitative variables. This step also produces the supplementary figure S2 file ('Figure_S2_pairs.pdf'). The second step consists in selecting the best P450 quantitative variables within the full model. Finally, a backward stepwise regression approach is used to select the model. 
+ **MyMetab_reg_NTSR.R:** this script described the regression of CYP6CY3 expression against metabolic resistant proportion. This consists in comparing different models with a maximum threshold value. This modelling is done only on homozygous sensitive clones. 
+ **MyMetab_test_TSRcomp.R:** this script is used to investigate the effect of the R81T genotype on the LC50 values, with or whithout PBO. 
+ **MyMetab_plot_Fig1_LD50.R:** the script to produce the Figure 1A 'correlation between LD50 with or without PBO' and Figure 1B 'LC50 dumbellplot with or without PBO'. 
+ **MyMetab_plot_Fig2_S1_copyexpr.R:** the script to produce the Figure 2 and Fugure S1 that illustrate the copy number and expression levels of CYP6 genes for each clones. 
+ **MyMetab_plot_Fig3_TSRcomp.R:** the script to produce the Figure 3 depicting the distribution of LC50 with PBO by TSR genotype
+ **MyMetab_plot_Fig4_NTSRprop.R:** the script to produce the Figure 4 depicting the evolution of the proportion of metabolic resistance as a function of CYP6CY3


## R session info
For reproducibility purpose, you will find all the information about the versions of R, Rstudio, OS etc., as well as the list and version number of the packages used at the time of publishing this script in the **session_info.txt** file.


## Citation
You can cite the related study as follow: 
+ Mottet M., Caddoux L., Fontaine S., Plantamp C., Bass C. and Barrès B. [*Myzus persicae* resistance to neonicotinoids - unravelling the contribution of different mechanisms to phenotype. *Pest Management Science*, 2024. doi.org/10.1002/ps.8316.](https://doi.org/10.1002/ps.8316)
 

If you want to use (some of) the code found on this page or if you want to cite this repository: 
+ Barrès B. [Supporting data and code for: *Myzus persicae* resistance to neonicotinoids - unravelling the contribution of different mechanisms to phenotype. Zenodo.](https://doi.org/10.5281/zenodo.10960421)