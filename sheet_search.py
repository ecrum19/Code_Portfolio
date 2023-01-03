import pandas as pd
import os
loc = open('/Users/eliascrum/PutontiLab/Ecoli_Proj/Ecolis.xlsx', 'r')

def excel_parse (location, column):
    df = pd.read_excel(location, sheet_name='patric3_query', usecols=column)  # reads xlsx
    output_list_x = df["Isolation Source"]      # assigns values of one column to list
    return output_list_x


def output(out, l):
    o = open(out, 'a')
    for i in l:
        o.write(str(i) + '\n')
    o.close()


output('Urine_Ecoli.txt', excel_parse(loc, 'AI'))
