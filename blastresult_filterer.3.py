#!/usr/bin/python3.6

import argparse, sys

parser=argparse.ArgumentParser()

parser.add_argument('infile', help='Input file name (in blastn outfmt 6 or 7)')
parser.add_argument('outfile', help='Name of output file to write')
parser.add_argument('--minidentity', help='Minimum accepted identity level of sequences (float)',
                    type = float, default = 50.0)
parser.add_argument('--maxevalue', help='Maximum accepted e-value (float)', type = float, default = 1e-10)
parser.add_argument('--minlen', help='Minimum alignment length (int). Takes priority over calcminlen',
                    type = int, default = 0)
parser.add_argument('--multi', help='Print multiple hits per query if true', action = 'store_true')
parser.add_argument('--calcminlen',
                    help='calculate minimum alignment length based on query length (query length * calcminlen) (float)',
                    type = float, default = 0)

args=parser.parse_args()

infile = args.infile
outfile = args.outfile
max_evalue = args.maxevalue
min_identity = args.minidentity
multi = args.multi
calcmin = args.calcminlen

def calc_min_len(queryID):
    """calculate minimum alignment length if this option is specified"""
    start = int(queryID.split('-')[1])*(-1)
    end = int(queryID.split(':')[1].split('-')[0])*(-1)
    minlen = (end - start) * calcmin
    return(minlen)

def get_minlen(argmin, argcalcmin, query):
    """Define minimum alignment length"""
    if argmin != 0:
        minlen = args.minlen
    elif argcalcmin != 0:
        minlen = calc_min_len(query)
    else:
        minlen = 200
    return(minlen)

prev_id = ""
prev_evalue = 1
queries = []
with open("blastresult_filterer.log", 'w') as out:
    out.write("FILTERED BLAST RESULTS\n")
    out.write(f"Printing queries and unique gene names with\n")
    out.write(f"\t - minimum identity of {min_identity}\n")
    out.write(f"\t - maximum e-value of {max_evalue}\n")
    if args.minlen != 0:
        out.write(f"\t - minimum alignment length: {args.minlen}")
    elif args.calcminlen != 0:
        out.write(f"\t - minimum alignment length calculated: query length * {args.calcminlen}")
    else:
        out.write(f"\t - minimum alignment length: defaulting to 200")

with open(outfile, 'w') as out:
    with open(infile) as f:
        for row in f:
            if row.startswith( "# Query" ) or "hits found" in row:
                out.write(row)
            elif not row.startswith("#"):
                array = row.split("\t")
                query = array[0]
                evalue = float(array[10])
                identity = float(array[2])
                alignlen = int(array[3])
                minlen = get_minlen(args.minlen, calcmin, query)
                if alignlen >= minlen and identity >= min_identity and evalue <= max_evalue:
                    if query not in queries:
                        queries.append(query)
                        first_evalue = evalue
                        out.write(f"{row.strip()}\n")
                    elif query in queries and multi == True and evalue <= first_evalue:
                        out.write(f"{row.strip()}\n")
