import os
intact_dir = '/Users/eliascrum/Programs/PythonProjects/Lab/Ecoli/intact_phage_seqs/'
new_intact_dir = '/Users/eliascrum/Programs/PythonProjects/Lab/Ecoli/new_intact_phage_seqs/'


def all_phages():
    def edit_distance(string1, string2):
        nope = False
        if len(string1) != len(string2):
            if len(string1) > len(string2) and len(string1) - len(string2) < max(0.01*len(string1), 0.01*len(string2)):
                difference = len(string1) - len(string2)
                string1 = string1[:len(string2)]
            elif len(string2) > len(string1) and len(string2) - len(string1) < max(0.01*len(string1), 0.01*len(string2)):
                difference = len(string2) - len(string1)
                string2 = string2[:len(string1)]
            else:
                nope = True
        else:
            difference = 0

        if not nope:
            for i in range(len(string1)):
                if string1[i] != string2[i]:
                    difference += 1
            return difference
        else:
            return 100000000000

    same = []
    for ar in os.listdir(intact_dir):
        if '.fasta' in ar:
            aa = open(intact_dir + ar, 'r')
            a = aa.read().split('\n')[1]

        for br in os.listdir(intact_dir):
            if '.fasta' in br:
                ba = open(intact_dir + br, 'r')
                b = ba.read().split('\n')[1]
                e = edit_distance(a, b)
                if e/((len(a) + len(b))/2) <= 0.01 and a != b:
                    same.append((ar, br))

        for cr in os.listdir(new_intact_dir):
            if '.fasta' in cr:
                ca = open(new_intact_dir + cr, 'r')
                c = ca.read().split('\n')[1]
                e = edit_distance(a, c)
                if e / ((len(a) + len(c)) / 2) <= 0.01 and a != c:
                    same.append((ar, cr))

    def new():
        for ds in os.listdir(new_intact_dir):
            if '.fasta' in ds:
                da = open(new_intact_dir + ds, 'r')
                d = da.read().split('\n')[1]

            for fs in os.listdir(intact_dir):
                if '.fasta' in fs:
                    fa = open(intact_dir + fs, 'r')
                    f = fa.read().split('\n')[1]
                    e = edit_distance(d, f)
                    if (e/len(d)) <= 0.01 and d != f:
                        same.append((ds, fs))

            for gs in os.listdir(new_intact_dir):
                if '.fasta' in gs:
                    ga = open(new_intact_dir + gs, 'r')
                    g = ga.read().split('\n')[1]
                    e = edit_distance(d, g)
                    if (e / len(d)) <= 0.01 and d != g:
                        same.append((ds, gs))

    o = open('/Users/eliascrum/Programs/PythonProjects/Lab/Ecoli/' + 'same_new.out', 'w')
    for a in same:
        o.write(str(a) + '\n')
    o.close()


def parse_outputs(filenames):
    path = '/Users/eliascrum/Programs/PythonProjects/Lab/Ecoli/'
    for f in filenames:
        a = open(path + f, 'r')
        ar = a.read().split('\n')
        all_same = {}
        for w in ar[:-1]:
            a_s = w[1:-1].split(',')
            if a_s[0] not in all_same:
                all_same[a_s[0]] = a_s[1]
            else:
                curr = all_same[a_s[0]]
                if type(curr) == list:
                    curr.append(a_s[1])
                else:
                    curr = [all_same[a_s[0]]]
                    curr.append(a_s[1])
                all_same[a_s[0]] = curr


    total = set()
    pointers = []
    for s in all_same:
        total.add(s)
        current = []
        current.append(s)
        if type(all_same[s]) == list:
            for k in all_same[s]:
                total.add(k)
                current.append(k)
        else:
            total.add(all_same[s])
            current.append(all_same[s])
        pointers.append(current)

    lens = []
    highest = 0
    h = []
    for g in pointers:
        lens.append(len(g))
        if len(g) > highest:
            highest=len(g)
            h = g
    print(sorted(lens)[-1])
    print(len(pointers))

    o = open('/Users/eliascrum/Programs/PythonProjects/Lab/Ecoli/' + 'same_all.out', 'w')
    for a in pointers:
        o.write(str(a) + '\n\n')
    o.close()


#all_phages()
parse_outputs(['same_new.out', 'same.out'])
