#!/usr/bin/env bash

"""
Author: Rauf Salamzade

Testing for developers of zol to make sure new code doesn't break anything
and to check proper functioning of the code.
"""

# Step 0: Uncompress test_case.tar.gz and cd into it.
rm -rf testing_cases/
tar -zxvf testing_cases.tar.gz
cd testing_cases/

############################################
### Aspergillus prepTG testing 
############################################

prepTG -i Aflavus_Genomes_GenBank/ -o prepTG_Aspergillus/

prepTG -i Aflavus_Genomes_FASTA/ -o prepTG_Aspergillus_Miniprot/ -rp GCF_009017415.1_ASM901741v1_protein.faa

##############################################
### Epa Enterococcus dataset based testing ###
##############################################

# Step 1a: run prepTG on target genomes directory (FASTA) 
prepTG -i Target_Genomes/ -o prepTG_Database/ -c 10 -cst -ma

# Step 2a: run fai (input type 1) to identify orthologous instances of epa from E. faecalis
fai -i Epa_MIBiG_GenBank.gbk Efaecium_EnGen0003_Epa_extraction.gbk -tg prepTG_Database/ -o fai_Results_1/ -c 10 --generate-plots -gdm HMM

# Step 2b: run fai (input type 2) to identify orthologous instances of epa from E. faecalis
fai -r Efaecalis_V583_Genome.fasta -rc NC_004668.1 -rs 2083902 -re 2115174 -tg prepTG_Database/ -o fai_Results_2/ -c 10 --generate-plots

# Step 2c: run fai (input type 3) to identify orthologous instances of epa from E. faecalis
fai -pq Epa_Proteins_from_MIBiG_GenBank.faa -tg prepTG_Database/ -o fai_Results_3/ -c 10 --generate-plots

# Step 3: run zol to perform comparative investigations of gene-clusters
ls -lht fai_Results_2/Final_Results/Homologous_Gene_Cluster_GenBanks/ | awk '{print $9}' | grep 'faecalis' > Efaecalis_GCs.txt 
zol -i fai_Results_2/Final_Results/Homologous_Gene_Cluster_GenBanks/ -o zol_Comparative_Results/ -c 10 --full-genbank-labels -f Efaecalis_GCs.txt
zol -i fai_Results_2/Final_Results/Homologous_Gene_Cluster_GenBanks/ -o zol_DC_Results/ -c 10 -dco -dcp "--approx-id 50 --mutual-cover 25"
zol -i fai_Results_2/Final_Results/Homologous_Gene_Cluster_GenBanks/ -o zol_InDepth_Results/ -c 10 -s -b -ces -aec
zol -i fai_Results_2/Final_Results/Homologous_Gene_Cluster_GenBanks/ -o zol_Derep_Reinflate_Results/ -c 10 -d -ri
zol -i fai_Results_2/Final_Results/Homologous_Gene_Cluster_GenBanks/ -o zol_Derep_IBC_Results/ -c 10 -d -ibc
zol -i fai_Results_2/Final_Results/Homologous_Gene_Cluster_GenBanks/ -o zol_Dom_Results/ -c 10 -dom
zol -i fai_Results_2/Final_Results/Homologous_Gene_Cluster_GenBanks/ -o zol_SFP_Results/ -c 10 -sfp
zol -i fai_Results_2/Final_Results/Homologous_Gene_Cluster_GenBanks/ -o zol_PFP_Custom_Results/ -c 10  -pfp="Full None 1e-5"

# Step 4: test out cgc and cgcg for collapsed gene cluster visualization generation
cgc -i zol_Comparative_Results/ -o cgc_Results/
cgcg -i zol_Comparative_Results/ -o cgcg_Results/

# Step 5: test out abon, atpoc, apos, salt, and zol-scape
abon -tg prepTG_Database/ -a antiSMASH_Results/Efaecalis_V583/ -g GECCO_Results/Efaecalis_V583/ -o abon_Results/ -c 10
atpoc -tg prepTG_Database/ -i Efaecalis_V583_Genome.fasta -ps PhiSpy_Results/Efaecalis_V583/ -vi VIBRANT_Results/Efaecalis_V583/ -gn geNomad_Results/Efaecalis_V583/ -o atpoc_Results/ -c 10
apos -tg prepTG_Database/ -i Efaecalis_V583_Genome.fasta -ms MOB-suite_Results/Efaecalis_V583/ -gn geNomad_Results/Efaecalis_V583/ -o apos_Results/ -c 10
salt -f fai_Results_1/ -tg prepTG_Database/ -o salt_Results/ -c 10
zol-scape -i Efaecalis_antiSMASH_Results/ -r Efaecalis_BiG-SCAPE_Results/ -o zol-scape_Results/ -c 10

# Step 6: test out regex for extracting genomic regions
mkdir -p regex_Results/

# Test 6a: Extract region from FASTA file (will run pyrodigal for gene calling)
regex -i Efaecalis_V583_Genome.fasta -s NC_004668.1 -b 2083902 -e 2115174 -o regex_Results/Efaecalis_epa_from_fasta.gbk

# Test 6b: Extract region from GenBank file (direct extraction)
regex -i Epa_MIBiG_GenBank.gbk -s BGC0000792 -b 1 -e 10000 -o regex_Results/Epa_subset_from_genbank.gbk

# Test 6c: Test with different gene calling methods on FASTA
regex -i Efaecalis_V583_Genome.fasta -s NC_004668.1 -b 2090000 -e 2100000 -o regex_Results/Efaecalis_region_prodigal.gbk -gcm prodigal

# Test 6d: Test with metagenomics mode
regex -i Efaecalis_V583_Genome.fasta -s NC_004668.1 -b 2095000 -e 2105000 -o regex_Results/Efaecalis_region_meta.gbk -m
