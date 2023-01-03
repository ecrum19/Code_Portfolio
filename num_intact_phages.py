import os


def new_phage_nums():
    path = '/media/catherine/ExtraDrive1/Eli_Ecoli/phaster_output/'
    all = os.listdir(path)

    list = open(path + 'num_intact.txt','w')
    no = ['num_intact.txt', 'all_questionable_phasters.fna', 'all_intact_phasters.fna', 'all_incomplete_phasters.fna']

    num = 0
    for i in all:
        if i not in no and '_output' in i:
            f = open(path+i+'/summary.txt')
            fa = f.readlines()
            for j in fa:
                if 'Totally' in j:
                    if '0' in j:
                        num += 1
                    list.write(i + ': ' + j + '\n')
    print(num)
    list.close()


def phage_nums_total(path, all):
    num = 0
    for i in all:
        ir = open(path+i, 'r')
        num += ir.read().count('>')
    print(num)

def len_phages(file):
    path = '/media/catherine/ExtraDrive1/Eli_Ecoli/phaster_output/'
    l = open(path + file, 'r')
    ls = l.read().strip().split('\n')
    num = 0
    for i in ls:
        if '>' in i:
            num += 1
    print(num)


def just_nums(f1):
    path = '/media/catherine/ExtraDrive1/Eli_Ecoli/phaster_output/'
    fi = open(path + f1, 'r')
    fis = fi.read().strip().split('\n')
    for d in fis:
        if d == '':
            fis.remove(d)

    of = open('intact_phages_guide.txt', 'w')

    all = {}
    for k in fis:
        name = k.split('_')[1]
        number = k.split(' ')[2]
        all[name] = number

    curr = 10
    while curr >= 0:
        for e in all:
            if int(all[e]) == curr:
                of.write('%s\t%s\n' % (e, all[e]))
        curr = curr - 1


def intact_phage_fastas(inf,p):
    f = open(inf, 'r')
    fo = f.read().strip().split('>')
    fo.pop(0)
    for s in fo:
        ss = s.split('\n')
        name = ss[0].replace('/', '')
        no = open(p+'/intact_phage_seqs/' + name + '.fasta', 'w')
        no.write('>' + s)
        no.close()


def phage_seq_lengths(inp, outp):
    out = open(outp+'/old_phagelengths.txt','w')
    lengths = {}
    jv = []
    for file in os.listdir(inp):
        o = open(inp+'/'+file, 'r')
        ore = o.read().split('>')
        ore.pop(0)
        ores = ore[0].split('\n')
        lengths[ores[0]] = len(ores[1])
        jv.append(len(ores[1]))
    for k in sorted(jv):
        for f in lengths:
            if lengths[f] == k:
                out.write('Length: %d\t%s\n' % (k, f))
    out.close()


def hist_phage_numbers():
    path = '/Users/eliascrum/Programs/PythonProjects/Lab/Ecoli/'
    old = open(path + 'num_intact.txt', 'r')
    new = open(path + 'Induction/num_intact.txt', 'r')

    o = old.read().split('\n\n')
    n = new.read().split('\n\n')

    intact_nums = open(path+'bact_intact_phage_counts.fna', 'w')

    for j in o[:-1]:
        intact_nums.write('%s\t%d\n' % (j.split(' ')[0][:-1], int(j.split(' ')[2])))
    for k in n[:-1]:
        intact_nums.write('%s\t%d\n' % (k.split(' ')[0][:-1], int(k.split(' ')[2])))
    intact_nums.close()


def big_multi_fasta():
    def remove_val(list, val):      # for removing values from splitting
        try:
            while True:
                list.remove(val)
        except ValueError:
            pass
        return list

    os.chdir('/Users/eliascrum/Programs/PythonProjects/Lab/Ecoli/')
    o = open('/Users/eliascrum/Programs/PythonProjects/Lab/Ecoli/all_prophage_genomes_short.fna.fna', 'r')

    w = open('new_prophage_multi.fna', 'w')
    g = remove_val(o.read().split('>'), '')
    for i in g:
        ih = i.split('\n')
        name = ''
        for h in ih[0].split('_')[1:]:
            if h != 'length':
                name += h + '_'
            else:
                break
        w.write('>' + name[:-1] + '\n' + ih[1] + '\n')

    w.close()


def old_multi_fasta():
    def remove_val(list, val):      # for removing values from splitting
        try:
            while True:
                list.remove(val)
        except ValueError:
            pass
        return list

    os.chdir('/Users/eliascrum/Programs/PythonProjects/Lab/Ecoli/')
    o = open('/Users/eliascrum/Programs/PythonProjects/Lab/Ecoli/all_prophage_genomes_short.fna', 'r')


    g = remove_val(o.read().split('>'), '')
    old = []
    for f in g:
        if 'AWS' not in f:
            old.append(f)

    c = 0
    while c < len(old):
        w = open('prophage_multi_%s.fna' % str(int(c/301)), 'w')
        for j in old[c:c+301]:
            js = j.split('\n')
            w.write('>' + js[0] + '\n' + js[1] + '\n')
        c += 301
        w.close()


def big_multi_fasta2():
    os.chdir('/Users/eliascrum/Programs/PythonProjects/Lab/Ecoli/')
    o = open('/Users/eliascrum/Programs/PythonProjects/Lab/Ecoli/prophage_multi_long.fna', 'r')
    b = open('/Users/eliascrum/Programs/PythonProjects/Lab/Ecoli/new_prophage_multi_short.fna', 'w')
    path='/Users/eliascrum/Programs/PythonProjects/Lab/Ecoli/new_intact_phage_seqs'
    for i in os.listdir(path):
        ih = open(path+'/'+i, 'r').read()
        o.write(ih)

    o.close()


def induced_phage_genomes_multifasta():
    files = os.listdir('/Users/eliascrum/Programs/PythonProjects/Lab/Ecoli/new_induced_phages/Induced_Phages')
    induced_phage_multi = open('/Users/eliascrum/Programs/PythonProjects/Lab/Ecoli/new_induced_phages/Induced_Phages/induced_phage_multi.fna', 'w')
    for f in files:
        if 'AWS' in f:
            o = open('/Users/eliascrum/Programs/PythonProjects/Lab/Ecoli/new_induced_phages/Induced_Phages/' + f, 'r')
            induced_phage_multi.write(o.read() + '\n')
    induced_phage_multi.close()


def induced_phage_genome_paths():
    files = os.listdir('/Users/eliascrum/Programs/PythonProjects/Lab/Ecoli/new_induced_phages/Induced_Phages')
    induced_phage_multi = open(
        '/Users/eliascrum/Programs/PythonProjects/Lab/Ecoli/new_induced_phages/Induced_Phages/induced_phage_paths.fna',
        'w')
    for f in files:
        if 'AWS' in f:
            induced_phage_multi.write('/Users/eliascrum/Programs/PythonProjects/Lab/Ecoli/new_induced_phages/Induced_Phages/' + f + '\n')


def induced_phage_blast():
    os.chdir('/media/catherine/ExtraDrive1/Eli_Ecoli/induced_phages_blast/')
    os.system('makeblastdb -in induced_phage_multi.fna -out induced -title induced -dbtype nucl')
    os.system('blastn -query %s -db induced -out induced_phages_out.csv -outfmt "10 qseqid sseqid pident length qstart qend sstart send evalue"'
        % ('/media/catherine/ExtraDrive1/Eli_Ecoli/induced_phages_blast/induced_phage_multi.fna'))



path = '/media/catherine/ExtraDrive1/Eli_Ecoli/new_phaster_output/'
all = ['all_new_intact_phasters.fna']

#phage_nums_total(path, all)

#induced_phage_genome_paths()
#induced_phage_blast()
#just_nums('num_intact.txt')
#len_phages('all_incomplete_phasters.fna')
path = '/Users/eliascrum/Programs/PythonProjects/Lab/Ecoli/Induction/'
#intact_phage_fastas(path + '/all_intact_phasters.fna', path)
#phage_seq_lengths(path+'/intact_phage_seqs', path)
hist_phage_numbers()
#old_multi_fasta()
#induced_phage_genomes_multifasta()

