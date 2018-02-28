#!/usr/bin/env python

"""
fastaParser.py

  Version 1.0 -- read in a multi-fasta file, parse, and output as required

"""

import sys, os, re, string
from optparse import OptionParser

usage = "usage: %prog [options]"
parser = OptionParser(usage=usage)
parser.add_option("-a", "--action", dest="action", 
                  help="Required query action, can be filter/adapter/convert")
parser.add_option("-i", "--input", dest="input", help="Input fasta filename, defaults to stdin.")
parser.add_option("-o", "--output", dest="output", help="Output filename, defaults to stdout.")
parser.add_option("-t", "--tdf", dest="tdf", default=False, action="store_true",
                  help="Option to set output file as tab delimited file, defaults to fasta format.")
parser.add_option("-w", "--width", dest="width", default=50, type='int',
                  help="Output fasta file width, defaults to 50bp each line.")
parser.add_option("-l", "--filter_len", dest="filter_len", default=100, type='int',
                  help="Filter sequences with length less then this value, defaults to 100bp.")
parser.add_option("-A", "--adapter", dest="adapter", help="Adapter sequence filename, fasta file, \
required for action adapter.")
parser.add_option("-s", "--adpt_min_len", dest="adpt_min_len", default=1, type='int',
                  help="Minimum length of the overlap between 3' end of a sequence and the \
adapter sequence, defaults to 1.")
parser.add_option("-v", "--verbose", dest="verbose", action="store_true", default=False)

(options, args) = parser.parse_args()

if not options.action:
    parser.print_help()
    sys.exit()


revcomp_table = string.maketrans("ACGT","TGCA")
def revcomp(seq):
    return string.translate(seq, revcomp_table)[::-1]


# a class for one sequence in a fasta file
class Sequence:
    def __init__(self):
        self.name = ""
        self.fields = []
        self.sequence = ""
    
    def assign(self, name_line):
        rec = name_line[1:].split(' ')
        self.name = rec[0]
        for i in range(1, len(rec)):
            if re.search(':', rec[i]):
                val = rec[i].split(':')
                self.fields.append(val[1])
            else:
                self.fields.append(rec[i])

    def append(self, seq):
        self.sequence += seq

    def length(self):
        return len(self.sequence)

    def empty(self):
        if len(self.sequence) == 0: return True
        return False

    def description(self):
        return string.join(self.fields,' ')

    def write(self, out_file):
        out_file.write(">%s %s\n" % (self.name,self.description()))
        width = options.width
        for i in range(0, self.length(), width):
            out_file.write("%s\n" % self.sequence[i:i+width])
        out_file.write("\n")

    def write_tdf(self,out_file):
        out_file.write("%s\t" % (self.name))
        for i in range(0, len(self.fields)):
            out_file.write("%s\t" % self.fields[i])
        out_file.write("%s\n" % self.sequence)

    def check_adapter(self, adapters):
        overlap_len = 0
        adp_name = ""
        for adpt in adapters.sequences:
            adp_seq = revcomp(adpt.sequence)
            min_len = options.adpt_min_len
            max_len = len(adp_seq)
            if self.length() < max_len: max_len = self.length()
            for overlap in range(max_len, min_len, -1):
                seq_suffix = self.sequence[-overlap:]
                adp_prefix = adp_seq[:overlap]
                if seq_suffix == adp_prefix:
                    space = ""
                    for i in range(0, self.length()-overlap): space += ' '
                    if overlap > overlap_len: 
                        overlap_len = overlap
                        adp_name = adpt.name
                    break
        if options.verbose and overlap_len > 0: 
            sys.stderr.write("%s\t%s\t%s\n" % (self.name, adp_name, overlap_len))
        return overlap_len

    def write_no_adapter(self, out_file, adp_len):
        width = options.width
        if self.length()-adp_len > 0:
            out_file.write(">%s %s\n" % (self.name, self.description()))
            if adp_len > 0:
                out_file.write("%s\n" % self.sequence[:-adp_len])
            else:
                out_file.write("%s\n" % self.sequence)


class Sequences:
    def __init__(self, fasta_file):
        self.sequences = []
        line = fasta_file.readline()
        while line != "":
            if line[0] == '>':
                seq = Sequence()
                seq.assign(line[:-1])
                self.sequences.append(seq)
            else:
                seq.append(line[:-1])
            line = fasta_file.readline()
    
    def write(self, out_file):
        for i in range(0, len(self.sequences)):
            seq = self.sequences[i]
            seq.write(out_file)


def fastaFilterByLength(input, output, min_len):
    seq = Sequence()
    line = input.readline()
    while line != "":
        if line[0] == '>':
            if seq.length() >= min_len:
                seq.write(output)
            seq = Sequence()
            seq.assign(line[:-1])
        else:
            seq.append(line[:-1])
        line = input.readline()
    if seq.length() >= min_len:
        seq.write(output)



def fastaFilterByAdapter(input, output, adapter_file):
    adp_file = open(adapter_file, 'r')
    adapters = Sequences(adp_file)
    adp_file.close()
    seq = Sequence()
    line = input.readline()
    while line != "":
        if line[0] == '>':
            if not seq.empty():
                adp_len = seq.check_adapter(adapters) 
                seq.write_no_adapter(output,adp_len)
            seq = Sequence()
            seq.assign(line[:-1])
        else:
            seq.append(line[:-1])
        line = input.readline()
    adp_len = seq.check_adapter(adapters)
    seq.write_no_adapter(output,adp_len)


def fastaConvert(input, output):
    seq = Sequence()
    line = input.readline()
    while line != "":
        if line[0] == '>':
            if not seq.empty():
                seq.write_tdf(output)
            seq = Sequence()
            seq.assign(line[:-1])
        else:
            seq.append(line[:-1])
        line = input.readline()
    seq.write_tdf(output)


###### Main #######
def main():
    if options.input: input = open(options.input, 'r')
    else: input = sys.stdin
    if options.output: output = open(options.output, 'w')
    else: output = sys.stdout
    if options.action == "filter":
        fastaFilterByLength(input,output,options.filter_len)
    elif options.action == "adapter":
        fastaFilterByAdapter(input,output,options.adapter)
    elif options.action == "convert":
        fastaConvert(input,output)
    if options.input:  input.close()
    if options.output: output.close()


main()

