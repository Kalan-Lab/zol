#!/usr/bin/env bash

#######################################################################################################
##### This is same test set as for conda based installation but using the run_ZOL.sh wrapper instead.
#######################################################################################################

# change permissions for run_ZOL.sh in this directory to allow execution.
chmod +x run_ZOL.sh

# Step 0: Download and uncompress test_case.tar.gz and cd into it.
wget https://github.com/Kalan-Lab/zol/blob/main/test_case.tar.gz
rm -rf test_case/
tar -zxvf test_case.tar.gz
cd test_case/

# Step 1: run prepTG on target genomes directory to prepare for fai
run_ZOL.py prepTG -i Target_Genomes/ -o prepTG_Database/ -c 4

# Step 2a: run fai (input type 1) to identify orthologous instances of epa from E. faecalis
run_ZOL.py fai -i Epa_MIBiG_GBK/ -tg prepTG_Database/ -o fai_Results_1/ -c 4 --generate_plots

# Step 2b: run fai (input type 2) to identify orthologous instances of epa from E. faecalis
run_ZOL.py fai -r Efaecalis_V583_Genome.fasta -rc NC_004668.1 -rs 2083902 -re 2115174 -tg prepTG_Database/ \
           -o fai_Results_2/ -c 4 --generate_plots

# Step 2c: run fai (input type 3) to identify orthologous instances of epa from E. faecalis
run_ZOL.py fai -pq Epa_Proteins_from_MIBiG_GenBank.faa -tg prepTG_Database/ -o fai_Results_3/ -c 4 --generate_plots

# Step 3: run zol to perform comparative investigations of gene-clusters
run_ZOL.py zol -i fai_Results_2/Homologous_GenBanks_Directory/ -o zol_Results/ -c 4 --full_genbank_labels
