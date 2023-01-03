import os
from Bio.Seq import Seq
import statistics


def remove_val(list, val):  # for removing values from splitting
    try:
        while True:
            list.remove(val)
    except ValueError:
        pass
    return list

def find_integrases():
    all = os.listdir('/Users/eliascrum/Programs/PythonProjects/Lab/Ecoli/Integrase')
    gb = []
    fna = []
    for i in all:
        if '.gb' in i:
            gb.append(i)
        elif '.fna' in i:
            fna.append(i)

    integrases_locations = {}
    for file in gb:
        cf = open('/Users/eliascrum/Programs/PythonProjects/Lab/Ecoli/Integrase/' + file, 'r')
        cgb = remove_val(cf.read().split('LOCUS'), '')
        for genome in cgb:
            if 'integrase' in genome or 'Integrase' in genome:
                info = genome.split('\n')
                name = remove_val(info[0].split(' '), '')[0]

                cds = genome.split('CDS')
                for s_cds in cds:
                    if 'integrase' in s_cds or 'Integrase' in s_cds:
                        l1 = s_cds.split('\n')[0]
                        location = remove_val(l1.split(' '), '')[-1]
                        if name in integrases_locations:
                            new_str = integrases_locations[name][0] + ',' + location
                            integrases_locations[name] = [new_str]
                        else:
                            integrases_locations[name] = [location]

    int_loc = open('/Users/eliascrum/Programs/PythonProjects/Lab/Ecoli/Integrase/integrase_locations.csv', 'w')
    for d in integrases_locations:
        int_loc.write(d + ',' + integrases_locations[d][0] + '\n')
    int_loc.close()


def get_integrase_seqs():
    locs = open('/Users/eliascrum/Programs/PythonProjects/Lab/Ecoli/Integrase/integrase_locations.csv', 'r')
    all_sequences = open('/Users/eliascrum/Programs/PythonProjects/Lab/Ecoli/all_prophage_genomes_short.fna', 'r')

    all_seqs = remove_val(all_sequences.read().split('>'), '')
    n_l = locs.read().split('\n')

    sequences = {}
    for integrases in n_l:
        name = integrases.split(',')[0]
        location = [i for i in integrases.split(',')[1:]]

        new_entry = False
        loc = 0
        while not new_entry and loc < len(all_seqs):
            if name in all_seqs[loc]:
                for num in location:
                    if 'complement' not in num:                 # to determine if location is on comp or not
                        start = num.split('..')[0]
                        end = num.split('..')[-1]
                        int_seq = all_seqs[loc].split('\n')[1][int(start):int(end) + 1]
                    else:
                        start = num.split('..')[0][11:]
                        end = num.split('..')[-1][:-1]
                        int_seq = str(Seq(all_seqs[loc].split('\n')[1][int(start):int(end) + 1]).reverse_complement())

                    if name not in sequences:                   # puts names and seqs into dict
                        sequences[name] = [int_seq]
                    else:
                        more = [g for g in sequences[name]]
                        more.append(int_seq)
                        sequences[name] = more
                    new_entry = True
            loc += 1

    output = open('/Users/eliascrum/Programs/PythonProjects/Lab/Ecoli/Integrase/integrase_sequences.fna', 'w')
    for seq in sequences:
        if len(sequences[seq]) > 1:
            for j in range(len(sequences[seq])):
                output.write('>' + seq + '_' + str(j) + '\n' + sequences[seq][j] + '\n')
        else:
            output.write('>' + seq + '\n' + sequences[seq][0] + '\n')
    output.close()


def cluster_info():
    os.chdir('/home/catherine/Desktop/')
    out = open('/home/catherine/Desktop/Eli/cluster_summary.csv', 'w')
    for file in os.listdir('/home/catherine/Desktop/Eli'):
        if 'fa' in file:
            f = open('/home/catherine/Desktop/Eli/' + file, 'r')
            fi = f.read().split('>')[1:]
            out.write('Cluster_%s,%d\n' % (file.split('.fa')[-1], len(fi)))
    out.close()


def integrase_lengths():
    l = open('/Users/eliascrum/Programs/PythonProjects/Lab/Ecoli/Integrase/good_integrase_seqs.fasta', 'r')
    o = open('/Users/eliascrum/Programs/PythonProjects/Lab/Ecoli/Integrase/good_integrase_sequence_lens.csv', 'w')
    nl = []
    for j in l.read().split('>')[1:]:
        name = j.split('\n')[0]
        seq = j.split('\n')[1]
        nl.append((len(seq), name))
    for k in sorted(nl):
        o.write('>' + k[1] + ',' + str(k[0]) + '\n')
    o.close()


def integrase_avg_len():
    l = open('/Users/eliascrum/Programs/PythonProjects/Lab/Ecoli/Integrase/bad_integrase_sequence_lens.csv', 'r')
    vals = []
    for j in l.read().split('\n')[:-1]:
        vals.append(int(j.split(',')[-1]))
    print(statistics.mean(vals), statistics.median(vals))


def filterIntegrases():
    new = open('/Users/eliascrum/Programs/PythonProjects/Lab/Ecoli/Integrase/shrunken_integrase_alignment.fasta', 'r')
    old = open('/Users/eliascrum/Programs/PythonProjects/Lab/Ecoli/Integrase/integrase_sequences.fasta', 'r')
    bad = open('/Users/eliascrum/Programs/PythonProjects/Lab/Ecoli/Integrase/bad_integrase_seqs.fasta', 'w')
    good = open('/Users/eliascrum/Programs/PythonProjects/Lab/Ecoli/Integrase/good_integrase_seqs.fasta', 'w')
    n = new.read().split('>')
    names = []
    for i in n[1:]:
        names.append(i.split('\n')[0])
    for j in old.read().split('>')[1:]:
        if j.split('\n')[0] not in names:
            bad.write('>' + j)
        else:
            good.write('>' + j)
    good.close()
    bad.close()


#filterIntegrases()
integrase_lengths()
integrase_avg_len()
#cluster_info()
#find_integrases()
#get_integrase_seqs()

# Cluster Command -- put usearch in path
# '/home/catherine/Software/usearch11.0.667_i86linux32' -cluster_fast /home/catherine/Desktop/Eli/integrase_sequences.fna -id 0.35 -centroids /home/catherine/Desktop/Eli/0.35_clusters/integrase_clusters.fasta

# kalign command:
# kalign infile outfile -format fasta

# FastTree Command:
# nohup FastTree -gtr -nt /Users/eliascrum/Programs/PythonProjects/Lab/Ecoli/Integrase/good_integrase_alignment.fasta > good_integrase_tree.out

