def rename():
    f = open('/Users/eliascrum/Programs/PythonProjects/Lab/Ecoli/new_prophage_taxonomy.csv', 'r')
    out = open('/Users/eliascrum/Programs/PythonProjects/Lab/Ecoli/all_coliphage_taxa.csv', 'w')

    fi = f.read().strip().split('\n')
    for i in fi:
        try:
            ni = i.split(',')
            if ni[1] in ['unclassified Caudovirales']:
                ni.pop(-1)
                ni.append('Unclassified Caudovirales')
            elif ni[1] not in ['Podoviridae', 'Siphoviridae', 'Myoviridae', 'Unclassified bacterial virus', ]:
                ni.pop(-1)
                ni.append('Unknown')
        except:
            continue
        name = ''
        for k in ni[0].split('_')[1:]:
            if k != 'length':
                name += k + '_'
            else:
                break
        out.write(name[:-1] + '\t' + ni[1] + '\n')
    out.close()


def combine_taxonomies():
    f = open('/Users/eliascrum/Programs/PythonProjects/Lab/Ecoli/all_coliphage_taxa.csv', 'r')
    f0 = open('/Users/eliascrum/Programs/PythonProjects/Lab/Ecoli/coliphage_taxa_full_corrected.txt', 'r')
    out = open('/Users/eliascrum/Programs/PythonProjects/Lab/Ecoli/total_coliphage_taxa.txt', 'w')

    fi = f.read().strip().split('\n')
    fi0 = f0.read().strip().split('\n')
    for i in fi:
        out.write(i + '\n')
    for j in fi0:
        out.write(j + '\n')

    out.close()

#rename()
combine_taxonomies()




