#!/usr/bin/python3.6

#Calculate FPKM from a table  with samples and features and a custom bed file
#table: column: sample
#       row: feature
#bed file: any feature region, has to include ID
#simplified calculation, for circRNA

import sys
import argparse
import pandas as pd

parser = argparse.ArgumentParser()
parser.add_argument('infile', help='tsv file of read counts')
parser.add_argument('features', help='feature file (bed)')
parser.add_argument('outfile')
args = parser.parse_args()  

circ = pd.read_csv(args.infile, sep = '\t', header = 0) #all read data in a csv file

with open(args.infile) as f:
    samples = f.readline().strip().split('\t')[1:] #sample name list

srsum = [ circ[s].sum() for s in samples] #list of sum total reads per sample

#feature length by bed file
flen = {}
with open(args.features) as ff:
    for line in ff:
        if not line.startswith("Region"):
            cid = line.strip().split('\t')[3]
            clen = abs(int(line.strip().split('\t')[2])-int(line.strip().split('\t')[1]))
            flen[cid] = int(clen)

with open(args.outfile, 'w') as out:
    with open(args.infile) as f:
        out.write(f.readline())
        for line in f:
            if not line.startswith('Region'):
                reads = line.strip().split('\t')[1:]
                cid = line.strip().split('\t')[0]
                out.write(f"{cid}")
                for i, read in enumerate(reads):
                    try:
                        #rpkm: (10^6 * C)/(N * L)
                        #in this case it's fpkm because my input counts supporting read pairs
                        fpkm = (1_000_000 * int(reads[i])) / (int(srsum[i]) * int(flen[cid]))
                    except:
                        #this happens only if there is no read in a sample
                        #prevent zero division error
                        fpkm = "NA"
                    out.write(f"\t{fpkm}")
                out.write(f"\n")

                
