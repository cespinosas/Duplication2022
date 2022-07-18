# Early effects of gene duplication on the robustness and phenotypic variability of gene regulatory networks

## Yuridia S. Posadas-García, Carlos Espinosa-Soto

Instituto de Fı́sica, Universidad Autónoma de San Luis Potosı́, Av. Parque Chapultepec 1570, 78295, San Luis Potosı́, Mexico Jul 18, 2022

## Description of contents

In these files you will find all the essential programs for reproducing the different experiments of this manuscript. All the programs are written in C++ and most of them require the GSL numeric library  (https://www.gnu.org/software/gsl/). The libraries and compilation flags required by `g++` are commented at the beginning of every program, as well as the necessary arguments. The programs are listed in the order in which they must be run to obtain the raw data.

### libs/ 

All the libraries with the functions needed for all the programs are in this directory. You need to compile these before anything else (`g++ -Wall -c *.cc`).

### outs/ 

This is the output directory. All the data from main.cc is stored inside this folders. Other programs look for the specific data to run in these files.

### bolsadebichos.cc 

This program obtains and stores the GRNs through a Monte Carlo walk.

### main.cc 

This program takes the GRNs and performs general gene addition and mutational experiments. Among other things, it collects the distance between phenotypes after gene addition and mutations. The program obtains similarity and mutational robustness.  This program also gets the number of accessible phenotypes. Other programs need the output files that this program produces.

### similarity_permut_pergene.cc 

This program performs experiments to assess the effect of different kinds of mutation (e.g. interaction-deleting mutations). 

### num_conex.cc 

This program counts the number of connections in the duplicate gene.

### dist_rec_phenotypes.cc 

This program detects and stores the recurrently accessible phenotypes.

### norec_phenotypes.cc 

This program store all the accessible phenotypes that are non-recurrent.

### spearmyphenotype.cc 

This program counts and classifies accessible phenotypes and the number of mutations leading to them. 

### mut_lead_rec.cc 

This program obtains the number of mutations leading to recurrently accessible phenotypes.

### sm_rec_norec_ori.cc 

This program calculates the similarity between the accessible phenotypes and the original phenotype.

### figures_and_ANOVA.R 

This script contains the essential code for the figures in the article and the type II ANOVA test. You will need to arrange the data as indicated.
