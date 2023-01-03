import os
import glob
import Bio
from Bio import SeqIO
import csv
import numpy as np

# Ec phage
list_seq = list(SeqIO.parse('/media/catherine/ExtraDrive1/Eli_Ecoli/pangenome_all.csv', "fasta"))

phage_dict = {}
for i in list_seq:
    line = i.id.split(("|"))
    gene_cluster = line[1][13:]
    genome_name = line[2][12:]
    print(genome_name)

    if genome_name not in phage_dict.keys():
        phage_dict[genome_name]=list()

    phage_dict[genome_name].append(gene_cluster)


genomes_list = list(phage_dict.keys())

values = list()
for i in range(len(genomes_list)):
    values.append(list())

outfile = open('edge_list_ecPhage_final.csv', 'w')
for i in range(len(genomes_list)):
    genome_x = genomes_list[i]
    genes_for_x = phage_dict[genome_x]
    for j in range(len(genomes_list)):
        genome_y = genomes_list[j]
        genes_for_y = phage_dict[genome_y]
        intersection = [value for value in genes_for_x if value in genes_for_y]
        values[i].append(len(intersection))

# write an edge list
for i in range(len(genomes_list)):
    for j in range(len(genomes_list)):
        if (int(str(values[i][j])) != 0):
            if (str(genomes_list[i] != str(genomes_list[j]))):
                # print(genomes_list[i])
                # print(","+ genomes_list[j])
                # print(','+str(values[i][j])+ "\n")
                outfile.write(genomes_list[i])
                outfile.write("," + genomes_list[j])
                outfile.write(',' + str(values[i][j]) + "\n")


# anvi-gen-genomes-storage -e Ec_phage_paths.txt -o EcPhage-GENOMES.db
# anvi-pan-genome -g EcPhage-GENOMES.db -n Ec_phage --mcl-inflation 2 --minbit 0.35 -T 25
# python3 phage_network.py

# cytoscape
