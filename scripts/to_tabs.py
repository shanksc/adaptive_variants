import sys
import pandas as pd
#convert space delimited file to tabs

csv = sys.argv[1]
out = sys.argv[2]

usage = '<csv> <out>'

df = pd.read_csv(csv)
#print(df.head())
#print(df.info())

contents = []
with open(csv, 'r') as f:
    for line in f.readlines():
        contents.append(line.strip().split())

with open(out, 'w') as f:
    for line in contents:
        f.write('\t'.join(line)+'\n')
