#!/usr/bin/env python

# New Tuxedo Workflow - Stage 1 of 3
# Align the RNA-seq reads to the genome

# Corresponds to steps 1 and 2 of following protocol:
# Pertea et al (2016) Nature Protocols, Vol. 11(9): 1650-1667; doi:10.1038/nprot.2016.095
# 1. Map the reads for each sample to the reference genome.
# 2. Sort and convert the SAM files to BAM.

# INPUT:
# idx       tar archive containing genome index files
# fastq     names (full path) of fastq files
# numcpus   number of CPUs available

# OUTPUT:
# BAM files             *.bam
# alignment statistics  *.alnstats

# Suggested settings for SBG platform (based on SBG HISAT2 tool)
# Number of CPUs: 8
# Memory: 31000MB

# python script to align RNA-seq reads
VERSION = "0.0.1"

import argparse
import re
import glob
from subprocess import call
    
def alignAssemble(args):
    
    call(["tar", "-xf", args.idx[0], "-C", "/idx"])
    idxPath = glob.glob('/idx/*.ht2')[0]
    idxPath = re.sub(r".\d.ht2$", "", idxPath)
    
    args.fastq.sort()
    
    reads1=[]
    for fq in args.fastq:
        match = re.search('.*_1.fastq.gz$', fq)
        if match:
            reads1 += [fq]

    reads2=[]
    for fq in args.fastq:
        match = re.search('.*_2.fastq.gz$', fq)
        if match:
            reads2 += [fq]
            
    samples=[]
    for s in reads1:
        s = re.sub(r"^/.*/", "", s)
        s = re.sub(r"_1.fastq.gz$", "", s)
        samples += [s]
    
    for i in range(0, len(samples), 1):
        cmd = "hisat2 -p " + str(args.numcpus[0]) + " --dta -x " + idxPath + " -1 " + reads1[i] + " -2 " + reads2[i] + " -S " + samples[i]+".sam 2>"+samples[i]+".alnstats"
        call(cmd, shell=True)
        cmd = "samtools sort -@ " + str(args.numcpus[0]) + " -o " + samples[i]+".bam " + samples[i]+".sam"
        call(cmd, shell=True)
        cmd = "rm " + samples[i] + ".sam"
       
if __name__ == "__main__":
    """ Parse the command line arguments """
    parser = argparse.ArgumentParser(description='Align RNA-seq reads.')
    parser.add_argument('idx', nargs=1, help='tar archive containing genome index files')
    parser.add_argument('fastq', metavar='FASTQ', nargs='+', help='names and paths of fastq files')
    parser.add_argument('--numcpus', type=int, metavar='NUMCPU', nargs=1, help='Number of CPUs available', default=[1])
    parser.add_argument("--version", action='version', version=VERSION) 
    args = parser.parse_args()

    """ Run the desired methods """
    alignAssemble(args)




