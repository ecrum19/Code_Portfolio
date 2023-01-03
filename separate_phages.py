import os
import glob
import Bio
from Bio import SeqIO
import csv
import numpy as np

'''
Workflow from phaster output file 'all_intact_phasters.fna' to a phage network.
'''


def sep_phasteroutput():
    f = open('/Users/eliascrum/Programs/PythonProjects/Lab/Ecoli/new_intact_phasters_1.fna', 'r')

    wanted = []
    w = open('/Users/eliascrum/Programs/PythonProjects/Lab/Ecoli/Potential_induced_phages.txt', 'r')
    wa = w.read().strip().split('\n')
    for a in wa:
        if 'trimmed' in a:
            wanted.append(a)
    phages = f.read().strip().split('>')
    try:
        phages.remove('')
    except:
        pass

    for i in phages:
        ni = i.split('\n')
        try:
            ni.remove('')
        except:
            pass
        if ni[0] in wanted:
            o = open('/Users/eliascrum/Programs/PythonProjects/Lab/Ecoli/new_induced_phages/' + ni[0][:36] + '.fa', 'w')
            o.write('>%s\n%s' % (ni[0], ni[1]))


def cleanupFastas():
    os.chdir('/media/catherine/ExtraDrive1/Eli_Ecoli/new_intact_phage_seqs')
    al = os.listdir('/media/catherine/ExtraDrive1/Eli_Ecoli/new_intact_phage_seqs')
    for k in al:
        os.system('anvi-script-reformat-fasta %s -o %s-f.fa -l 0 --simplify-names --report-file problems.txt' % (k, k[:-3]))
        os.system('rm %s' % k)


def makedb():
    os.chdir('/media/catherine/ExtraDrive1/Eli_Ecoli/new_intact_phage_seqs')
    ai = os.listdir('/media/catherine/ExtraDrive1/Eli_Ecoli/new_intact_phage_seqs')
    for f in ai:
        os.system("anvi-gen-contigs-database -f %s -o /media/catherine/ExtraDrive1/Eli_Ecoli/new_coliphage_db/%s.db -n 'Coliphage_contig'" % (f, f[:-3]))


def hmmAndCogs():
    os.chdir('/media/catherine/ExtraDrive1/Eli_Ecoli/new_coliphage_db')
    dl = os.listdir('/media/catherine/ExtraDrive1/Eli_Ecoli/new_coliphage_db')
    for d in dl:
        os.system('anvi-run-hmms -c %s --num-threads 12' % d)
        os.system('anvi-run-ncbi-cogs -c %s --num-threads 12' % d)


def cleannames():
    os.chdir('/media/catherine/ExtraDrive1/Eli_Ecoli/all_coliphage_db')
    ni = os.listdir()

    path = '/media/catherine/ExtraDrive1/Eli_Ecoli/all_coliphage_db/'
    for w in ni:
        if w[-2:] == '-f':
            os.system('mv %s%s %s%s' % (path, w, path, w[:-2]))


def one_dir():
    curr_path = '/media/catherine/ExtraDrive1/Eli_Ecoli/new_coliphage_db'
    os.chdir(curr_path)
    new_path = '/media/catherine/ExtraDrive1/Eli_Ecoli/all_coliphage_db/'
    for w in os.listdir(curr_path):
        os.system('mv %s/%s %s%s' % (curr_path, w, new_path, w))

    curr_path = '/media/catherine/ExtraDrive1/Eli_Ecoli/coliphage_db'
    os.chdir(curr_path)
    for w in os.listdir(curr_path):
        os.system('mv %s/%s %s%s' % (curr_path, w, new_path, w))


def makeAnvioGenomelist():
    paths = os.listdir('/media/catherine/ExtraDrive1/Eli_Ecoli/all_coliphage_db')

    pathlist = open('/media/catherine/ExtraDrive1/Eli_Ecoli/all_Ec_phage_paths.txt', 'w')
    pathlist.write('name\tcontigs_db_path\n')
    for j in paths:
        pathlist.write(j + '\t')
        pathlist.write('/media/catherine/ExtraDrive1/Eli_Ecoli/all_coliphage_db/' + j + '\n')
    pathlist.close()


def phage_network():            # Need to have populated Anvi'o Pan-Genome file
    # coliphage
    list_seq = list(SeqIO.parse('/media/catherine/ExtraDrive1/Eli_Ecoli/all_pangenome_all.csv', "fasta"))

    phage_dict = {}
    for i in list_seq:
        line = i.id.split(("|"))
        gene_cluster = line[1][13:]
        genome_name = line[2][12:]
        print(genome_name)

        if genome_name not in phage_dict.keys():
            phage_dict[genome_name] = list()

        phage_dict[genome_name].append(gene_cluster)

    genomes_list = list(phage_dict.keys())

    values = list()
    for i in range(len(genomes_list)):
        values.append(list())

    outfile = open('edge_list_final.csv', 'w')
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




#sep_phasteroutput()    # parses all intact phages found into separate files
#cleanupFastas()        # removes characters Anvio doesn't like via anvi-script-reformat-fasta
#makedb()               # makes .db files from cleaned fasta files
#hmmAndCogs()           # annotates dbs
#cleannames()           # for anvi-gen-genomes-storage (it doesn't like .db)
#one_dir()
#makeAnvioGenomelist()  # to call as an argument for anvi-gen-genomes-storage

os.chdir('/media/catherine/ExtraDrive1/Eli_Ecoli')

os.system('anvi-gen-genomes-storage -e all_Ec_phage_paths.txt -o all_EcPhage-GENOMES.db')
os.system('anvi-pan-genome -g all_EcPhage-GENOMES.db -n all_Ec_phage --mcl-inflation 2 --minbit 0.35 -T 25 --skip-hierarchical-clustering')
os.system('anvi-get-sequences-for-gene-clusters -p all_Ec_phage/all_Ec_phage-PAN.db -g all_EcPhage-GENOMES.db -o all_pangenome_all.csv')

#phage_network()

# or if anvi'o can't make a pan-genome (because there are no core genes or too many genomes)

# conda activate anvio-6.2 (or whatever version is installed)
# Rscript --vanilla /media/catherine/ExtraDrive1/Eli_Ecoli/netprep.r -d /media/catherine/ExtraDrive1/Eli_Ecoli/all_Ec_phage
# (Ec_phage is folder created from anvi-pan-genome command)

'''
Anvio Workflow:
0. Separate phage genomes + clean them (remove non ASCII characters)
1. anvi-gen-contigs-database
2. anvi-gen-genomes-storage
3. anvi-pan-genome
4. phage network script
5. Cytoscape to visualize network
'''
