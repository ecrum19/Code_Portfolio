import os


def mke_dir_lst(d):         # makes list of files from one directory
    files = os.listdir(d)
    dirs = open('new_ecoli_genomes_lst.txt', 'w')
    for f in files:
        dirs.write(f + '\n')
    dirs.close()


def reform(d):
    g = open(d, 'a')
    for f in files:
        print(f + '\t' + '/data/Eli/GG_db/' + f)
    g.close()
    
    
def rename(d):       # renames anvio database files to exclude '.db' (because anvio doesn't like it for anvi-pan command
    files = os.listdir(d)
    for i in files:         # Renames files in dir with first 13 chars
        if len(i) > 13:
            os.rename(r'/Users/eliascrum/fastANI/Lacto_G/' + i, r'/Users/eliascrum/fastANI/Lacto_G/' + i[0:13] + '.fna')


def reciprocal():           # checks if there are reciprocal edges in a phage network edge file
    f = open('/Users/eliascrum/Programs/PythonProjects/Lab/Ecoli/genome_edgelist1.txt', 'r')
    edges = f.read().strip().split('\n')
    print('Edges: ' + str(len(edges)))
    nodes = {}
    edge_pairs = []
    for i in edges:
        curr_edge = i.split('\t')
        edge_pairs.append(curr_edge)
        if curr_edge[0] not in nodes:
            nodes[curr_edge[0]] = 1
        else:
            nodes[curr_edge[0]] = nodes[curr_edge[0]] + 1

        if curr_edge[1] not in nodes:
            nodes[curr_edge[1]] = 1
        else:
            nodes[curr_edge[1]] = nodes[curr_edge[1]] + 1

    print('Nodes: ' + str(len(nodes)))
    val = edge_pairs[125][1] + '\t' + edge_pairs[125][0]
    if val in edges:
        print('Duplicate: ' + val)
    return nodes


def numbers():  # prints number of lines in a file
    h = open('/media/catherine/ExtraDrive1/Eli_Ecoli/eli_ecoli_genomes_lst.txt','r')
    hs = h.readlines()
    print(len(hs))


def taxa_IDs(nodes):
    # takes in the output of a local BLAST (containing taxonomic info)
    # outputs a list of genome names and their taxonomic info (including 'Unknown' if no BLAST result)

    h = open('/Users/eliascrum/Programs/PythonProjects/Lab/Ecoli/Eli_coliphages_taxa.csv', 'r')
    out = open('/Users/eliascrum/Programs/PythonProjects/Lab/Ecoli/coliphage_taxa_full.txt', 'w')

    hs = h.read().strip().split('\n\n')

    have_taxa = []
    taxa = {}
    for q in nodes:
        for ha in hs:
            has = ha.split(',')
            if q in ha:
                have_taxa.append(q)
                taxa[q] = has[1]
                if 'Unclassified ' in taxa[q]:
                    taxa[q] = 'Unknown'
            elif q not in taxa:
                taxa[q] = 'Unknown'

    for a in taxa:
        out.write(a + '\t' + taxa[a] + '\n')

    out.close()

    # input the taxa doc to Cytoscape (use shared name)


# Driver
taxa_IDs(reciprocal())

#mke_dir_lst('/media/catherine/ExtraDrive1/Eli_Ecoli/ec_pt1')
