[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.10960421.svg)](https://doi.org/10.5281/zenodo.10960421)

# Supporting data and code for: *Myzus persicae* resistance to neonicotinoids - unravelling the contribution of different mechanisms to phenotype.
*This repository contains the R code used for the data analyses and production of the figures of the related article*

![alt text](https://db3pap005files.storage.live.com/y4mDlIdRciVaiIk5fs2Dl1bMvduiUuhuSt5MaPMCfoKs6iG_Wv1pzwiD5bzqd8Cc0khCc3R2vKX_Dy0G7SA_ZAKBkK7Z6ErM82rmN0PBrGVhi5J_X1MgkOaWnPbpG-9XsdJlQBP_EH7NkxoujphX2IPsmmx_-znhZq-RUVToP9UaiIKUqIxeI-bNseRTtwMZXvy?width=1584&height=588&cropmode=none)




## Context
Insecticide resistance can lead to the repeated application of treatments that are both ineffective in pest control and harmful for the environment. In order to establish sound strategies to slow down or limit these resistance phenomena, precise knowledge of the underlying mechanisms, their relative contribution and their potential interaction is essential. *Myzus persicae* is a major pest aphid, capable of infesting a wide range of crops, with significant economic impact. *M. persicae* is known to have evolved resistances to insecticides of different families, including the abundantly used neonicotinoids. *M. persicae* resistance to neonicotinoids has previously been described as being due to two main mechanisms: a P450 overexpression metabolic resistance and a target-site mutation, R81T, but their respective contribution to resistant phenotypes remains unclear. In this study we combined extensive bioassays with synergist on numerous clones, with gene copy number and expression quantification of two key P450 enzymes (*CYP6CY3* and *CYP6CY4*) to explain the observed phenotypes and assess the relative contribution of metabolic and target-site mechanisms to neonicotinoid resistance in *M. persicae*.


## Datasets
There are 3 datasets used in this study. The files can be found in the "data" folder. 

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
  + *clone-ID*: 
  + *genetic-group*: 
  + *nAChR-81*: 
  + *LC50*: 
  + *LC50-PBO*: 
  + *CY3_CN*: 
  + *CY4_CN*: 
  + *CY23_CN*: 
  + *CY3_EXP*: 
  + *CY4_EXP*: 
  + *CY23_EXP*: 
  + *CY3_SE_CN*: 
  + *CY4_SE_CN*: 
  + *CY23_SE_CN*: 
  + *CY3_SE_EXP*: 
  + *CY4_SE_EXP*: 
  + *CY23_SE_EXP*: 

+ **summaDataResc.txt:** the same dataset as '*summaData.txt*' but rescaled on another clone
  + *bla1*: fsdfsdfgsdfgf
  
  
## R scripts
+ **MyMetab_load.R:** the script to load the different datasets and packages in the environment
+ **MyMetab_cibleR_reg.R:** 
+ **MyMetab_drc.R:** 
+ **MyMetab_glm.R:**
+ **MyMetab_metabR_reg.R:** 
+ **MyMetab_plot_glm.R:** 
+ **MyMetab_plot_rawdat.R:** 

## Citation
You will soon be able (hopefully) to cite the related study as follow: 

