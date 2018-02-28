#!/usr/bin/env python


import sys, os, re, time
from optparse import OptionParser

from yzUtil import system_run
from yzSample import Sample
from yzBASE import *
from yzSAM import *

usage = "usage: %prog [options]"
parser = OptionParser(usage=usage)
parser.add_option("-i","--input",dest="input",
                  help="Input is a file containing a list of bam files.")
parser.add_option("-b","--bam",dest="bam",
                  help="Direct input bam files from command line, comma delimited. Exclusive with -i.")
parser.add_option("-f","--refasta",dest="refasta",
                  help="Required input reference fasta filename, need to have fai index file with it.")
parser.add_option("-o","--output",dest="output",
                  help="Option to set output filename, defaults to stdout.")
parser.add_option("-r","--region",dest="region",
                  help="Option to call SNPs in this region only.")
parser.add_option("-u","--qual",dest="qual",type='int',default=30,
                  help="Option to set the minimum mapping quality for SNP call. Defaults to 30.")
parser.add_option("-V","--vcf",dest="vcf",
                  default=False, action="store_true",
                  help="Option to output vcf files. Defaults to call/coverage/purity format.")
parser.add_option("-I","--indel",dest="indel",
                  default=False, action="store_true",
                  help="Option to call INDELs too.")
parser.add_option("-H","--hets",dest="hets",
                  default=False, action="store_true",
                  help="Option to call heterozygous SNPs too.")
parser.add_option("-S","--pos",dest="pos_file",
                  help="Option to provide the SNP positions from a file to get genotype values. \
The file needs to be tab delimited with two columns: seqid and coordinates.")
parser.add_option("-a","--param",dest="param",
                  help="Option to set samtools mpileup parameters.")
parser.add_option("-c","--min_cov",dest="cov",
                  default=3, type='int',
                  help="Option to set the minimum coverage for snp call, defaults to 3.")
parser.add_option("-p","--min_purity",dest="purity",
                  default=0.98, type='float',
                  help="Option to set the minimum purity for snp call, defaults to 0.98.")
parser.add_option("-d","--hets_cov",dest="hets_cov",
                  default=10, type='int',
                  help="Option to set the minimum coverage for hets call, defaults to 10.")
parser.add_option("-q","--hets_purity",dest="hets_purity",
                  default=0.7, type='float',
                  help="Option to set the maximum purity for hets call, defaults to 0.7.")
parser.add_option("-C","--indel_cov",dest="indel_cov",
                  default=1, type='int',
                  help="Option to set the minimum coverage for indel call, defaults to 1.")
parser.add_option("-P","--indel_purity",dest="indel_purity",
                  default=0.5, type='float',
                  help="Option to set the minimum purity for indel call, defaults to 0.5.")
parser.add_option("-Q","--indel_perc",dest="indel_perc",
                  default=0.1, type='float',
                  help="Option to set the minimum percent for displaying an indel allele, \
defaults to 0.1.")
parser.add_option("-s","--bsa",dest="bsa",
                  default=False, action="store_true",
                  help="Option to print Varients for Bulk Segregation Analysis (BSA) with \
the purity given by option -q.")
parser.add_option("-v","--verbose",dest="verbose",
                  default=False, action="store_true")

(options, args) = parser.parse_args()
if not options.refasta or (options.input and options.bam):
    parser.print_help()
    sys.exit(1)


def get_samples():
    samples = []
    if options.bam: 
        rec = options.bam.split(',')
        for i in range(0, len(rec)):
#            name = "sample%s" % (i+1)
            filename = rec[i]
            val = filename.split('/')
            name = '.'.join((val[-1].split('.'))[:-1])
            if os.path.exists(filename):
                samples.append(Sample(name,filename,'BAM'))
            else:
                sys.stderr.write("[bam2vcf.py] WARNING: skip %s, because it does not exist.\n" % filename)
    else:
        if options.input: in_file = open(options.input, 'r')
        else: in_file = sys.stdin
        line = in_file.readline()
        while line != "":
            if line[0] != '#':
                rec = line[:-1].split('\t')
                filename = rec[1]
                if os.path.exists(filename):
                    samples.append(Sample(rec[0],filename,'BAM'))
                else:
                    sys.stderr.write("[bam2vcf.py] WARNING: skip %s, because it does not exist.\n" % filename)
            line = in_file.readline()
        if options.input: in_file.close()
    return samples


def header_write(outfile, samples):
    outfile.write("#SeqID\tPosition\tType\tRefbase")
    for i in range(0, len(samples)):
        name = samples[i].name
        outfile.write("\t%s_Basecall\t%s_Cov\t%s_Purity" % (name,name,name))
        if options.bsa: outfile.write("\t%s_ACGT" % name)
    outfile.write("\n")


def bam_snpcall(samples, out_file):
    if out_file: outfile = open(out_file, 'w')
    else: outfile = sys.stdout
    header_write(outfile, samples)
    cmd = "samtools mpileup -f %s -q %s" % (options.refasta,options.qual)
    if options.region: cmd += " -r %s" % options.region
    if options.param: cmd += " %s" % options.param
    for i in range(0, len(samples)):
        cmd += " %s" % samples[i].file
    if options.verbose: sys.stderr.write("%s\n" % cmd)
    infile = os.popen("%s" % cmd)
    line = infile.readline()
    while line != "":
        if line[0] != '#':
            basepair = SAMpileupLineParser(line)
            if options.bsa and basepair.is_varient(options.cov, options.hets_purity): 
                basepair.write_varient(outfile)
            elif basepair.is_snp(options.cov, options.purity):
                basepair.write_snp(outfile,options.hets_cov,options.hets_purity)
            elif options.hets and basepair.is_hets(options.hets_cov, options.hets_purity):
                basepair.write_hets(outfile,options.hets_cov,options.hets_purity)
            if options.indel and basepair.is_indel(options.indel_cov, options.indel_purity):
                basepair.write_indel(outfile,options.indel_perc)
        line = infile.readline()
    if out_file: outfile.close()
    return
    

def sample_bcfcall(bamfile, vcffile):
    cmd = "samtools mpileup -uf %s" % options.refasta
    if options.region: cmd += " -r %s" % options.region
    if options.param:  cmd += " %s" % options.param
    cmd += " %s | bcftools call -v -m " % (bamfile)
    cmd += " - -o %s" % vcffile
    system_run(cmd)


def union_snps(vcfs, outfile):
    tmp_file = "%s.tmp" % outfile
    for i in range(0, len(vcfs)):
        cmd = "cat %s | grep -v \"#\"" % vcfs[i]
        if not options.hets: cmd += " | grep \"1/1\""
        cmd += " >> %s" % tmp_file
        system_run(cmd)
    cmd = "cat %s | cut -f1,2 | sort | uniq > %s" % (tmp_file, outfile)
    system_run(cmd)
    os.remove(tmp_file)
    

def GetSNPpositions(pos_file):
    positions = {}
    pos_fp = open(pos_file, 'rb')
    pos_line = pos_fp.readline()
    while pos_line != '':
        if pos_line[0] != '#':
            rec = pos_line[:-1].split('\t')
            seqid = rec[0]
            coord = int(rec[1])
            if not seqid in positions: positions[seqid] = []
            positions[seqid] += [coord]
        pos_line = pos_fp.readline()
    pos_fp.close()
    for seqid in positions:
        positions[seqid].sort()
    return positions


def GetRefbases(pos_file, refasta):
    refbases = {}
    cmd = '/mnt/fasta_extract %s -p -i %s -l 0 | /mnt/fastaParser.py -a convert' % (refasta, pos_file)
    if options.verbose: sys.stderr.write(cmd + '\n')
    fp = os.popen(cmd)
    line = fp.readline()
    while line != '':
        rec = line[:-1].split('\t')
        pos = rec[0]
        refbase = rec[1]
        refbases[pos] = refbase
        line = fp.readline()
    fp.close()
    return refbases

        
def PrintEmptyLine(outfile, seqid, coord, num_samples, refbases, delimit):
    outfile.write('%s%s%s%s%s' % (seqid, delimit, coord, delimit, 'SNP'))
    pos = '%s:%s' % (seqid,coord)
    if pos in refbases: refbase = refbases[pos]
    else: refbase = 'N'
    outfile.write('%s%s' % (delimit, refbase))
    for i in xrange(num_samples):
        outfile.write('%s%s' % (delimit, '.'))
        outfile.write('%s%s' % (delimit, 0))
        outfile.write('%s%.1f' % (delimit, 0.0))
    outfile.write('\n')


def snarp(samples, pos_file, outfile):
    if options.vcf: cmd = "samtools mpileup -uf %s" % options.refasta
    else: cmd = "samtools mpileup -f %s" % options.refasta
    if options.region: cmd += " -r %s" % options.region
    if options.qual:   cmd += " -q %s" % options.qual
    if options.param:  cmd += " %s" % options.param
    cmd += " -l %s" % pos_file
    for i in range(0, len(samples)):
        cmd += " %s" % samples[i].file
    if options.vcf: 
        cmd += " | bcftools view -"
        if outfile: cmd += " > %s" % outfile
        system_run(cmd)
    else:
        # Get SNP positions and refbases
        positions = GetSNPpositions(pos_file)
        refbases  = GetRefbases(pos_file, options.refasta)
        # Get coverage data from mpileup
        if options.verbose: sys.stderr.write("%s\n" % cmd)
        if outfile: out_fpt = open(outfile, 'wb')
        else: out_fpt = sys.stdout
        header_write(out_fpt, samples)
        visited = {}
        for seq_id in positions: visited[seq_id] = False
        fpt = os.popen(cmd)
        pos_seq_id = ''
        pos_line_id = 0
        line = fpt.readline()
        while line != '':
            if line[0] != '#':
                basepair = SAMpileupLineParser(line)
                if basepair.seqid != pos_seq_id:
                    if pos_seq_id:
                        while pos_line_id < len(positions[pos_seq_id]):
                            PrintEmptyLine(out_fpt, pos_seq_id, positions[pos_seq_id][pos_line_id], len(samples), refbases, '\t')
                            pos_line_id += 1
                    pos_seq_id = basepair.seqid
                    pos_line_id = 0
                    visited[pos_seq_id] = True
                while basepair.pos > positions[basepair.seqid][pos_line_id]:
                    PrintEmptyLine(out_fpt, basepair.seqid, positions[basepair.seqid][pos_line_id], len(samples), refbases, '\t')
                    pos_line_id += 1
                if basepair.pos == positions[basepair.seqid][pos_line_id]:
                    pos_line_id += 1
                if options.bsa: basepair.write_varient(out_fpt)
                else: basepair.write(out_fpt, options.hets_cov, options.hets_purity)
            line = fpt.readline()
        if pos_seq_id != '':
            while pos_line_id < len(positions[pos_seq_id]):
                PrintEmptyLine(out_fpt, pos_seq_id, positions[pos_seq_id][pos_line_id], len(samples), refbases, '\t')
                pos_line_id += 1
        for seq_id in visited:
            if visited[seq_id] == False:
                for i in xrange(len(positions[seq_id])):
                    PrintEmptyLine(out_fpt, seq_id, positions[seq_id][i], len(samples), refbases, '\t')
                visited[seq_id] = True
        if outfile: out_fpt.close()
    return


def vcf_snpcall(samples, outfile):
    vcfs = []
    for i in range(0, len(samples)):
        sample = samples[i]
        vcfile_name = "%s.vcf" % sample.name
        vcfs.append(vcfile_name)
        sample_bcfcall(sample.file, vcfile_name)
    if outfile: pos_file = "%s.pos" % outfile
    else: pos_file = "%s.pos" % int(time.time())
    union_snps(vcfs, pos_file)
    snarp(samples, pos_file, outfile)
    os.remove(pos_file)
    for filename in vcfs: os.remove(filename)


def main():
    samples = get_samples()
    if options.pos_file:
        snarp(samples, options.pos_file, options.output)
    elif options.vcf:
        vcf_snpcall(samples,options.output)
    else:
        bam_snpcall(samples,options.output)


if __name__ == '__main__':
    main()

