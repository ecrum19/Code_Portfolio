import os
path = '/media/catherine/ExtraDrive1/Eli_Ecoli/ecoli_genomes/superinfection_blast'


def fix_fastas():
    for a in os.listdir('/media/catherine/ExtraDrive1/Eli_Ecoli/ecoli_genomes/bridget_ecoli'):
        oa = open('/media/catherine/ExtraDrive1/Eli_Ecoli/ecoli_genomes/bridget_ecoli/' + a, 'r')
        ov = oa.read().strip().split('\n')
        nice_out = open('/media/catherine/ExtraDrive1/Eli_Ecoli/ecoli_genomes/bridget_ecoli/' + 'better_' + a, 'w')

        nice_g = ''
        for b in ov:
            if '>' in b:
                curr_g = '\n\n' + b + '\n'
                nice_g = nice_g + curr_g
            else:
                nice_g = nice_g + b
        nice_out.write(nice_g)
        nice_out.close()


def make_database():
    out = open('blast_ecs.fasta', 'w')
    for i in os.listdir('/media/catherine/ExtraDrive1/Eli_Ecoli/ecoli_genomes/superinfection_blast'):
        o = open(path + '/' + i, 'r')
        val = o.read().strip().split('\n')
        for j in val:
            if '>' in j or j == '':
                continue
            else:
                out.write('>%s\n%s\n' % (i[0:20], j))   # figure this out

    out.close()
    os.system('makeblastdb -in blast_ecs.fasta -out supinfecdb -title supinfecdb -dbtype nucl')


def run_blast():
    for j in os.listdir('/media/catherine/ExtraDrive1/Eli_Ecoli/supinfec_blast/Phages_of_interest'):
        os.system('blastn -query /media/catherine/ExtraDrive1/Eli_Ecoli/supinfec_blast/Phages_of_interest/%s -db supinfecdb -out %s_output.csv -outfmt "6 qseqid sacc pident length qstart qend sstart send evalue bitscore"' % (j, j))


def pull_results():
    l = ['AWS667_S32_L002_contigs_1_NO.fa_output.csv', 'AWS700_S39_L002_contigs_4_NO.fa_output.csv', 'AWS723_S48_L002_contigs_4_NO.fa_output.csv', 'AWS775_S61_L002_contigs_2_NO.fa_output.csv','AWS777_S62_L002_contigs_2_NO.fa_output.csv', 'AWS831_S75_L002_contigs_2_NO.fa_output.csv', 'AWS831_S75_L002_contigs_9_NO.fa_output.csv']
    m = ['i527.fa_output.csv', 'i6653.fa_output.csv', 'i6721.fa_output.csv']

    for file in m:
        os.system('scp catherine@10.23.19.202:/media/catherine/ExtraDrive1/Eli_Ecoli/supinfec_blast/%s /Users/eliascrum/Programs/PythonProjects/Lab/Ecoli/SuperInfectionCheck' % file)


#fix_fastas()
#make_database()
#run_blast()
pull_results()  # run on local terminal -- bioinfo16


#makeblastdb -in /media/catherine/ExtraDrive1/Eli_Ecoli/supinfec_blast/Phages_of_interest/i6653.fa -out i6653db -title i6653db -dbtype nucl
#blastn -query /media/catherine/ExtraDrive1/Eli_Ecoli/supinfec_blast/Phages_of_interest/i6721.fa -db i6653db -out %s_output.csv -outfmt "6 qseqid sacc pident length qstart qend sstart send evalue bitscore"
