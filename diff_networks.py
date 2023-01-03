import os

def network_sets(edge_list):
    edges = open('/Users/eliascrum/Programs/PythonProjects/Lab/Ecoli/network_2/' + edge_list, 'r')
    ae = edges.read().strip().split('\n')
    print(len(ae))
    networks = []
    a, b, c, d, e, f, g, h, i, k = set(), set(), set(), set(), set(), set(), set(), set(), set(), set()
    sets = [a, b, c, d, e, f, g, h, i, k]

    # initialize first set
    starta = ae[0].split('\t')
    a.add(starta[0])
    a.add(starta[1])

    remaining_nodes = ae
    iterations = 0

    def fill_new_set(s, t):
        for j in range(len(sets)):
            if len(sets[j]) == 0:
                sets[j].add(s)
                sets[j].add(t)
                break
            else:
                print('set %s contains nodes' % j)

    # first round
    while len(remaining_nodes) > 0 and iterations < 5:
        not_in = []
        pos = 0
        while pos < len(remaining_nodes):
            si = remaining_nodes[pos].split('\t')
            source = si[0]
            target = si[1]
            if source in a or target in a:
                a.add(target)
                a.add(source)
                pos += 1
            else:
                not_in.append(remaining_nodes[pos])
                pos += 1
        remaining_nodes = not_in
        iterations += 1

    # After First Round
    while len(remaining_nodes) > 0:
        not_in = []
        pos = 0
        while pos < len(remaining_nodes):
            si = remaining_nodes[pos].split('\t')
            source = si[0]
            target = si[1]
            for i in sets:
                if source in i or target in i:
                    i.add(target)
                    i.add(source)
                    pos += 1
                    break
            else:
                not_in.append(remaining_nodes[pos])
                pos += 1

        if len(remaining_nodes) == len(not_in):
            si = not_in.pop(0).split('\t')
            source = si[0]
            target = si[1]
            fill_new_set(source, target)

        remaining_nodes = not_in
        iterations += 1

    while len(sets) > 0:
        lar = 0
        curr_high = 0
        for a in range(len(sets)):
            curr_high = max(curr_high, len(sets[a]))
            if curr_high == len(sets[a]):
                lar = a
        networks.append(sets[lar])
        sets.pop(lar)

    return networks


def out(n):
    os.chdir('/Users/eliascrum/Programs/PythonProjects/Lab/Ecoli/network_2/')
    labels = ['a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j']
    o = open('all_network_guide.txt', 'w')
    for d in range(len(n)):
        o.write('Component %s\n' % labels[d])
        o.write(str(n[d]) + '\n' + str(len(n[d])))
        o.write('\n\n')
    o.close()


def nodes_left_out(edge_list, genomes):
    path = '/Users/eliascrum/Programs/PythonProjects/Lab/Ecoli/network_2/'
    edges = open(path + edge_list, 'r')
    names = open(path + genomes, 'r')
    in_network = set()
    for l in edges.read().split('\n')[:-1]:
        ls = l.split('\t')
        one = ls[0]
        two = ls[1]
        in_network.add(one)
        in_network.add(two)

    tot_network = set()
    for n in names.read().split('\n'):
        ns = n.split('\t')[0]
        tot_network.add(ns)

    for o in tot_network:
        if o not in in_network:
            print(o)

    return None


#out(network_sets('all_genome_edgelist3_weight.txt'))
nodes_left_out('all_genome_edgelist3_weight.txt', 'total_coliphage_taxa.txt')