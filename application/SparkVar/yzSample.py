#!/opt/common/compbio/bin/python

"""
yzSample.py
  A module for sample information.
"""

import re, os, sys


FileTypes = ['BAM','CDB','FILTAG','FASTA','FASTQ']

class Sample:
    def __init__(self,name,file,type):
        self.name = name
        self.file  = file
        self.type = type

    def write(self,outfile):
        outfile.write("%s\t%s\t" % (self.name,self.file))
        outfile.write("%s\n" % self.type)


def getSamples(infile, type='BAM'):
    samples = []
    line = infile.readline()
    line_id = 1
    while line != "":
        if line[0] != '#':
            rec = line[:-1].split('\t')
            if len(rec) < 2: 
                sys.stderr.write("[yzSample.py] ERROR: at least two columns provided.\n%s : %s" % (line_id, line))
                sys.exit(1)
            name = rec[0]
            file = rec[1]
            if len(rec) > 2: type = rec[2]
            samples.append(Sample(name,file,type))
        line_id += 1
        line = infile.readline()
    return samples


def writeSamples(samples, outfile):
    outfile.write("###### %s samples ######\n" % len(samples))
    for i in range(0, len(samples)):
        sample = samples[i]
        outfile.write("%s\t%s\t%s\n" % (sample.name,sample.file,sample.type))
    outfile.write("\n")


