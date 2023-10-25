import sys
#convert space delimited file to tabs

csv = sys.argv[1]
out = sys.argv[2]

usage = '<csv> <out>'

contents = []
with open(csv, 'r') as f:
    for line in f.readlines():
        contents.append(line.strip().split())

with open(out, 'w') as f:
    for line in contents:
        f.write('\t'.join(line)+'\n')
