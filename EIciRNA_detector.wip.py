#!/usr/bin/python3.6

#Goal: check intron retention
#Method:
#   Find circRNA specific reads for given circRNA from all samples
#   Check average intronic and exonic coverage
#      separately for all included introns and exons
#To do:
#   check why it didn't find find some reads that were supporting according to circstar
#   optimize stuff
#        - more exact read position limits
#   set up limits for intron coverage accepted as retention


import sys
import argparse
import pysam
import re
import subprocess

#Arguments
parser=argparse.ArgumentParser()
parser.add_argument('infile',
                    help = 'Combined bamfile of circRNA read trios')
parser.add_argument('outfile',
                    help = 'Output bam suffix for circRNA specific bam files. \
                    Prefix will be the circRNA ID (region)')
parser.add_argument('--region',
                    help = 'Specify genomic region for one circRNA: chr:start-end, one per line \
                    Only use if interested in only one region',
                    type = str)
parser.add_argument('--regions',
                    help = 'filename of circRNA region list: chr:start-end, one per line.')
parser.add_argument('--introns',
                    help = 'filename of introns list: chr:start-end, one per line.')
parser.add_argument('--exons',
                    help = 'filename of exon list: chr:start-end, one per line.')
parser.add_argument('--remove_bam',
                    help = 'if true, removes per circRNA generated bam files',
                    action = 'store_true')
args = parser.parse_args()

#Combined bam file of every circRNA read trio as input
f = pysam.AlignmentFile(args.infile, 'rb')

#Making exon and intron lists
with open(args.introns) as fi:
    ilist = [ line.rstrip('\n') for line in fi ]
with open(args.exons) as fe:
    elist = [ line.strip('\n') for line in fe ]

#Functions to use
def get_supporting_reads(reg):
   """Find trios with all 3 parts within region"""
   qdict = {}
   region = region_setup(reg)
   for aread in f.fetch(region = region):
      qname = aread.query_name
      if qname not in qdict.keys():
         new = {qname: 1}
         qdict.update(new)
      else:
         n = int(qdict[qname]) + 1
         qdict[qname] = n
   return(qdict)

def check_actual_support(reg, qdict):
   """Making sure read trios are actually supporting"""
   #basic setup
   supp = []
   region = region_setup(reg)
   start = region.split(":")[1].split("-")[0]
   end = region.split("-")[1]
   #check actual support
   for aread in f.fetch(region = region):
      qname = aread.query_name
      try:
         if qdict[qname] >= 3: #check if all reads are in region
            if aread.reference_start <= int(start) + 3 \
            and aread.reference_start >= int(start) - 3 \
            and re.search('[S,H]', aread.cigarstring):
               qdict[qname] += 1
            elif aread.reference_end <= int(end) + 2 \
            and aread.reference_end >= int(end) - 2 \
            and re.search('[S,H]', aread.cigarstring):
               qdict[qname] += 1
      except:
         pass
      #check if read is supporting
      if qdict[qname] == 5:
         supp.append(qname)
   return(supp)

def to_print(supp, circ):
   """Select what to output to circRNA specific bam"""
   for aread in f.fetch(region = region_setup(reg)):
      if aread.query_name in supp:
         circ.write(aread)

def region_setup(region):
   """check and fix occasional position switchups"""
   chrom = region.split(":")[0]
   pos1 = int(region.split(":")[1].split("-")[0])
   pos2 = int(region.split("-")[1])
   if pos1 < pos2:
      return(region)
   elif pos1 > pos2:
      region = f"{chrom}:{pos2}-{pos1}"
      return(region)   

def region_split(region):
   """check and fix occasional position switchups"""
   #Somewhat redundant, will have to combine with region_setup
   chrom = region.split(":")[0]
   pos1 = int(region.split(":")[1].split("-")[0])
   pos2 = int(region.split("-")[1])
   if pos1 < pos2:
       start = pos1
       end = pos2
   else:
       start = pos2
       end = pos1
   circreg = [chrom, start, end]
   return(circreg)

def get_introns(circreg, ilist):
   """find intronic or exonic regions in circRNA"""
   #Works for both, named after introns because i first did it for that
   inc_int = []
   for intron in ilist:
       ichrom = intron.split(":")[0]
       istart = int(intron.split(":")[1].split("-")[0])
       iend = int(intron.split("-")[1])
       if ichrom == circreg[0]:
           if istart >= circreg[1] - 1 \
           and iend <= circreg[2]:
               inc_int.append(intron)
   return(inc_int)


def run_samtools_depth(region):
   """Running samtools depth for given region"""
   covlist = []
   tmp = open("tmp.txt", 'w')
   subprocess.run(['samtools', 'index', cirbam])
   subprocess.run(['samtools', 'depth', '-aa', '-r', region, cirbam], stdout = tmp) #bookmark-comment

def calc_avg_cov():
   """Calculcating average depth for region"""
   with open("tmp.txt", 'r') as tmp:
       covlist = [ int(line.rstrip().split('\t')[2]) for line in tmp ]
   if len(covlist) != 0:
       avg = sum(covlist) / len(covlist)
       return(avg)

def make_covlist(introns):
   """Listing depth results"""
   icovs = []
   for i in introns:
       run_samtools_depth(i)
       icovs.append(calc_avg_cov())
   return(icovs)

def print_this(icovs, ecovs):
   """What to print as final result"""
   try:
       ecov_avg = sum(ecovs)/len(ecovs)
   except:
       ecov_avg = 0
   if ecov_avg != 0:
        #       print(ecov_avg)
        #       print(icovs)
       for i in icovs:
           if i / ecov_avg >= 0.4:
               out = f"{reg}\t{icovs}\t{ecovs}\t{ecov_avg}\tgood"
           elif i / ecov_avg >= 0.1:
               out = f"{reg}\t{icovs}\t{ecovs}\t{ecov_avg}\tmid"
           else:
               out = f"{reg}\t{icovs}\t{ecovs}\t{ecov_avg}\tlow"
   else:
       out = f"{reg}\t{icovs}\t0\t0\tnon-EI"
   try:
       return(out)
   except:
       print(f"error at {reg}")

#Run read selection main
if args.region: #For only one region
   reg = args.region
   circ = pysam.AlignmentFile(args.outfile, 'wb', template = f)
   qdict = get_supporting_reads(reg)
   supp = check_actual_support(reg, qdict)
   to_print(supp, circ)
   circ.close()

   cirbam = args.outfile
   circreg = region_split(reg)
   introns = get_introns(circreg, ilist)
   exons = get_introns(circreg, elist)
   icovs = make_covlist(introns)
   ecovs = make_covlist(exons)
   print(print_this(icovs, ecovs))  
  

elif args.regions: #For multiple regions
   with open(args.regions) as rf:
      regs = [ line.strip() for line in rf ]

   for reg in regs:
      #Start with full bam
      r = reg.strip().split(":")
      fname = f"{r[0]}-{r[1]}."
      circ = pysam.AlignmentFile(fname+args.outfile, 'wb', template = f)
      qdict = get_supporting_reads(reg)
      supp = check_actual_support(reg, qdict)
      to_print(supp, circ)
      circ.close()
      
      #continue working with circRNAspecific bam
      cirbam = fname+args.outfile
      circreg = region_split(reg)
      introns = get_introns(circreg, ilist)
      exons = get_introns(circreg, elist)
      icovs = make_covlist(introns)
      ecovs = make_covlist(exons)
      outp = print_this(icovs, ecovs)
      print(outp)
   
      if args.remove_bam:
          subprocess.run(['rm', cirbam])
          try:
              subprocess.run(['rm', cirbam+'.bai'])
          except:
              pass

f.close()

