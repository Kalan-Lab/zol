#!/usr/bin/env bash

# Step 0: Uncompress test_case.tar.gz and cd into it.
rm -rf test_case/
tar -zxvf test_case.tar.gz 
cd test_case/

# Step 1: run prepTG on target genomes directory to prepare for fai
zol-suite prepTG -i Target_Genomes/ -o prepTG_Database/ -c 4 -cst

# Step 2a: run fai (input type 1) to identify orthologous instances of epa from E. faecalis
zol-suite fai -i Epa_MIBiG_GBK/Epa_MIBiG_GenBank.gbk -tg prepTG_Database/ -o fai_Results_1/ -c 4 --generate-plots -gdm HMM

# Step 2b: run fai (input type 2) to identify orthologous instances of epa from E. faecalis
zol-suite fai -r Efaecalis_V583_Genome.fasta -rc NC_004668.1 -rs 2083902 -re 2115174 -tg prepTG_Database/ -o fai_Results_2/ -c 4 --generate-plots

# Step 2c: run fai (input type 3) to identify orthologous instances of epa from E. faecalis
zol-suite fai -pq Epa_Proteins_from_MIBiG_GenBank.faa -tg prepTG_Database/ -o fai_Results_3/ -c 4 --generate-plots

# Step 3: run zol to perform comparative investigations of gene-clusters
zol-suite zol -i fai_Results_2/Final_Results/Homologous_Gene_Cluster_GenBanks/ -o zol_Results/ -c 4 --full-genbank-labels

# Step 4: test out cgc and cgcg for collapsed gene cluster visualization generation
zol-suite cgc -i zol_Results/ -o cgc_Results/

zol-suite cgcg -i zol_Results/ -o cgcg_Results/
