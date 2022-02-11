import sys
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('infile', help ='single line peptid fasta file')

args = parser.parse_args()

circid = args.infile.split('.')[0].replace("-", ":", 1)
#print(circid)

#calculate ORF lengths
orf_lens = []
with open(args.infile) as inf:
    for line in inf:
        if not line.startswith(">"):
            orf_lens.append(len(line.strip()))

#Find maximum ORF length and all sequences this length
maxlen_nt = int(max(orf_lens)) * 3
maxlen_seqs = [i for i, e in enumerate(orf_lens) if e == max(orf_lens)]

#Gather peptid sequences and IDs
pepseq_ids = []
pepseqs = []
for seq in maxlen_seqs:
    pepseq_ids.append((2 * int(seq)) - 1)
    pepseqs.append(2 * int(seq))

#Screen input file and print important info
pepsids = []
peps = []
with open(args.infile) as inf:
    for i, line in enumerate(inf):
        if i in pepseq_ids:
            pepsids.append(line.strip())
        if i in pepseqs:
            peps.append(line.strip())

for a, z in enumerate(pepsids):
    print(f"{circid}\t{max(orf_lens)}\t{maxlen_nt}\t{pepsids[a]}\t{peps[a]}")
