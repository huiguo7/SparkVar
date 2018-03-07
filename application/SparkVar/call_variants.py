#!/usr/bin/env python


import sys
import os
import re
import time
from datetime import datetime
from optparse import OptionParser, OptionGroup
from BasePair import BasePair

PROGRAM_NAME = 'call_variants'


### system call ###
def syscall(cmd, verbose):
    if verbose: sys.stderr.write('(%s) [%s] #cmd: %s\n' % (str(datetime.now()), PROGRAM_NAME, cmd))
    retval = os.system(cmd)
    if retval:
        if verbose: sys.stderr.write('[%s] ERROR: system call failed: %s\n' % (PROGRAM_NAME, cmd))
        sys.exit(1)


### parse args and options ###
def parse_args():
    defaults = {'map_qual':30, 'min_cov':3, 'min_purity':0.98, \
                'hets_cov':10, 'hets_purity':0.7, 'indel_cov':1, 'indel_purity':0.5, 'indel_perc':0.1}
    usage = "usage: %prog [options]"
    parser = OptionParser(usage=usage)
    parser.add_option("-b","--bam_files",dest="bam_files",
                  help="Option to take input bam files from command line, comma delimited. Exclusive with -B.")
    parser.add_option("-B","--bam_list",dest="bam_list",
                  help="Option to take the inputs from a file of bam filename list, exclusive with -b.")
    parser.add_option("-f","--refasta",dest="refasta",
                  help="Required input reference fasta filename, need to have fai index file with it.")
    parser.add_option("-o","--output",dest="output",
                  help="Option to set output filename, defaults to stdout.")
    parser.add_option("-v","--verbose",dest="verbose",default=False, action="store_true")

    bam = OptionGroup(parser, 'Samtools mpileup Options')
    bam.add_option("-r","--region",dest="region",
                  help="Option to call variants in this region only.")
    bam.add_option("-u","--qual",dest="qual",type='int',default=defaults['map_qual'],
                  help="Option to set the minimum mapping quality, defaults to %s." % defaults['map_qual'])
    bam.add_option("","--pos",dest="pos_file",
                  help="Option to provide the variant positions from a file to get genotype values. \
The file needs to be tab delimited with two columns: seqid and coordinates.")
    bam.add_option("","--param",dest="param",
                  help="Option to set additional samtools mpileup parameters.")
    parser.add_option_group(bam)

    method = OptionGroup(parser, 'Variant Call Method Options')
    method.add_option("","--indel",dest="call_indel",default=False, action="store_true",
                  help="Option to call INDELs too.")
    method.add_option("","--hets",dest="call_hets",default=False, action="store_true",
                  help="Option to call heterozygous SNPs too.")
    method.add_option("","--vcf",dest="vcf",default=False,action="store_true",
                  help="Option to use bcftools to call variants, defaults to call by samtools mpileup.")
    parser.add_option_group(method)

    cutoff = OptionGroup(parser, 'Variant Call Cutoff Options')
    cutoff.add_option("","--min_cov",dest="cov",default=defaults['min_cov'],type='int',
                  help="Option to set the minimum coverage for snp call, defaults to %s." % defaults['min_cov'])
    cutoff.add_option("","--min_purity",dest="purity",default=defaults['min_purity'],type='float',
                  help="Option to set the minimum purity for snp call, defaults to %s." % defaults['min_purity'])
    cutoff.add_option("","--hets_cov",dest="hets_cov",default=defaults['hets_cov'],type='int',
                  help="Option to set the minimum coverage for hets call, defaults to %s." % defaults['hets_cov'])
    cutoff.add_option("","--hets_purity",dest="hets_purity",default=defaults['hets_purity'],type='float',
                  help="Option to set the maximum purity for hets call, defaults to %s." % defaults['hets_purity'])
    cutoff.add_option("","--indel_cov",dest="indel_cov",default=defaults['indel_cov'],type='int',
                  help="Option to set the minimum coverage for indel call, defaults to %s." % defaults['indel_cov'])
    cutoff.add_option("","--indel_purity",dest="indel_purity",default=defaults['indel_purity'],type='float',
                  help="Option to set the minimum purity for indel call, defaults to %s." % defaults['indel_purity'])
    cutoff.add_option("","--indel_perc",dest="indel_perc",default=defaults['indel_perc'],type='float',
                  help="Option to set the minimum percent for displaying an indel allele, defaults to %s." % defaults['indel_perc'])
    parser.add_option_group(cutoff)


    (options, args) = parser.parse_args()
    if not options.refasta or (options.bam_list and options.bam):
        parser.print_help()
        sys.exit(1)

    return options, args



#####################################
## Get Sample BAM files
#####################################
def get_samples(bamfiles, bam_list, delimit='\t'):
    samples = {}
    if bamfiles: 
        for filename in bamfiles.split(','):
            name = os.path.splitext(os.path.basename(filename))[0]
            if os.path.exists(filename):
                samples[name] = filename
            else:
                sys.stderr.write("[%s] WARNING: skip %s, because it does not exist.\n" % (PROGRAM_NAME,filename))
    else:
        if bam_list: in_file = open(bam_list, 'r')
        else: in_file = sys.stdin
        for line in in_file.readlines():
            if line[0] != '#':
                rec = line[:-1].split(delimit)
                (name, filename) = rec[0], rec[1]
                if os.path.exists(filename):
                    samples[name] = filename
                else:
                    sys.stderr.write("[%s] WARNING: skip %s, because it does not exist.\n" % (PROGRAM_NAME,filename))
        if bam_list: in_file.close()
    return samples



## Print Header
def HeaderLine(names, delimit):
    fields = ['CHR','COORDINATE','TYPE','REFBASE']
    for i in range(0, len(names)):
        name = names[i]
        fields += ['%s_BASECALL'%name,'%s_COV'%name,'%s_PURITY'%name]
    return delimit.join(fields)+'\n'


######################################
# variants call by samtools mpileup
######################################

class Variant_Cutoffs:
    def __init__(self):
        self.min_cov = 3           # minimum coverage for homo SNP call
        self.min_purity = 0.98     # minimum purity for homo SNP call
        self.hets_cov = 10         # minimum coverage for hets SNP call
        self.hets_purity = 0.7     # minimum purity for hets SNP call
        self.indel_cov = 1         # minimum coverage for INDEL call
        self.indel_purity = 0.5    # minimum purity for INDEL call
        self.indel_perc = 0.1      # minimum percent for displaying an INDEL allele

    def load_from_options(self, options):
        self.min_cov = options.cov
        self.min_purity = options.purity
        self.hets_cov = options.hets_cov
        self.hets_purity = options.hets_purity
        self.indel_cov = options.indel_cov
        self.indel_purity = options.indel_purity
        self.indel_perc = options.indel_perc


def variant_call_by_mpileup(samples, mpileup_cmd, out_file, cutoffs, call_hets, call_indel, delimit, verbose):
    names = samples.keys()
    names.sort()

    if out_file: outfile = open(out_file, 'w')
    else: outfile = sys.stdout
    outfile.write(HeaderLine(names, delimit))
    cmd = mpileup_cmd
    for i in range(0, len(names)):
        cmd += " %s" % samples[names[i]]
    if verbose: sys.stderr.write("%s\n" % cmd)
    infile = os.popen("%s" % cmd)
    line = infile.readline()
    while line != "":
        if line[0] != '#':
            rec = line[:-1].split(delimit)
            basepair = BasePair(rec)
            if basepair.is_snp(cutoffs.min_cov, cutoffs.min_purity):
                basepair.write_snp(outfile, cutoffs.hets_cov, cutoffs.hets_purity)
            elif call_hets and basepair.is_hets(cutoffs.hets_cov, cutoffs.hets_purity):
                basepair.write_hets(outfile, cutoffs.hets_cov, cutoffs.hets_purity)
            if call_indel and basepair.is_indel(cutoffs.indel_cov, cutoffs.indel_purity):
                basepair.write_indel(outfile, cutoffs.indel_perc)
        line = infile.readline()
    if out_file: outfile.close()
    return
    


########################################
#### Use mpileup to do snarping
########################################
def GetSNPpositions(pos_file, delimit='\t'):
    positions = {}
    pos_fp = open(pos_file, 'rb')
    pos_line = pos_fp.readline()
    while pos_line != '':
        if pos_line[0] != '#':
            rec = pos_line[:-1].split(delimit)
            seqid = rec[0]
            coord = int(rec[1])
            if not seqid in positions: positions[seqid] = []
            positions[seqid] += [coord]
        pos_line = pos_fp.readline()
    pos_fp.close()
    for seqid in positions:
        positions[seqid].sort()
    return positions


def GetRefbases(pos_file, refasta, verbose):
    refbases = {}
    cmd = 'fasta_extract %s -p -i %s -l 0' % (refasta, pos_file)
    if verbose: sys.stderr.write(cmd + '\n')
    fp = os.popen(cmd)
    line = fp.readline()
    snp_id = None
    while line != '':
        if line[0] == '>':
            snp_id = line[1:-1].strip()
            line = fp.readline()
            refbases[snp_id] = line[0]
        line = fp.readline()
    fp.close()
    return refbases

        
def EmptyLine(seqid, coord, num_samples, refbases, delimit):
    fields = [seqid, str(coord), 'UNK']
    pos = '%s:%s' % (seqid,coord)
    if pos in refbases: refbase = refbases[pos]
    else: refbase = 'N'
    fields += [refbase]
    for i in range(0, num_samples):
        fields += ['.','0','0.0']
    return delimit.join(fields)+'\n'


def snarp(samples, refasta, pos_file, mpileup_cmd, out_filename, delimit, verbose):
    # Get SNP positions and refbases
    positions = GetSNPpositions(pos_file)
    refbases  = GetRefbases(pos_file, refasta, verbose)

    # Initialize
    visited = {}
    for seq_id in positions: visited[seq_id] = False
    pos_seq_id, pos_line_id = '', 0
    names = samples.keys()
    names.sort()

    # open output file
    if out_filename: outfile = open(out_filename, 'w')
    else: outfile = sys.stdout
    outfile.write(HeaderLine(names, delimit))

    # call samtools mpileup
    cmd = mpileup_cmd + " -l %s" % pos_file
    for i in range(0, len(names)):
        cmd += " %s" % samples[names[i]]
    infile = os.popen(cmd)
    if verbose: sys.stderr.write("[%s] #cmd: %s\n" % (PROGRAM_NAME,cmd))

    # parse mpileup lines to call variants
    line = infile.readline()
    while line != '':
        if line[0] != '#':
            basepair = BasePair(line[:-1].split(delimit))
            if basepair.seqid != pos_seq_id:
                if pos_seq_id:
                    while pos_line_id < len(positions[pos_seq_id]):
                        outfile.write(EmptyLine(pos_seq_id, positions[pos_seq_id][pos_line_id], len(names), refbases, delimit))
                        pos_line_id += 1
                pos_seq_id = basepair.seqid
                pos_line_id = 0
                visited[pos_seq_id] = True
            while basepair.pos > positions[basepair.seqid][pos_line_id]:
                outfile.write(EmptyLine(basepair.seqid, positions[basepair.seqid][pos_line_id], len(names), refbases, delimit))
                pos_line_id += 1
            if basepair.pos == positions[basepair.seqid][pos_line_id]:
                pos_line_id += 1
            basepair.write(outfile, delimit)
        line = infile.readline()

    # print for variants with no read coverage
    if pos_seq_id != '':
        while pos_line_id < len(positions[pos_seq_id]):
            outfile.write(EmptyLine(pos_seq_id, positions[pos_seq_id][pos_line_id], len(names), refbases, delimit))
            pos_line_id += 1
    for seq_id in visited:
        if visited[seq_id] == False:
            for i in range(0, len(positions[seq_id])):
                outfile.write(EmptyLine(seq_id, positions[seq_id][i], len(names), refbases, delimit))
            visited[seq_id] = True

    if out_filename: outfile.close()
    return



########################################
#### Use bcftools to call variants  ####
########################################

# get the union of variant coordinates list
def vcf_union_snps(vcfs, outfile, call_hets, verbose):
    tmp_file = "%s.tmp" % outfile
    for i in range(0, len(vcfs)):
        cmd = "cat %s | grep -v \"#\"" % vcfs[i]
        if not call_hets: cmd += " | grep \"1/1\""
        cmd += " >> %s" % tmp_file
        syscall(cmd, verbose)
    cmd = "cat %s | cut -f1,2 | sort | uniq > %s" % (tmp_file, outfile)
    syscall(cmd, verbose)
    os.remove(tmp_file)
    

# call variants by bcftool
def variant_call_by_bcftools(samples, mpileup_cmd, out_file, call_hets, delimit, verbose):
    vcfs = []
    names = samples.keys()
    names.sort()
    for i in range(0, len(names)):
        name = names[i]
        vcf_name = "%s.vcf" % name
        vcfs.append(vcf_name)
        cmd = mpileup_cmd + ' %s | bcftools call -v -m - -o %s' % (samples[name],vcf_name)
        syscall(cmd, verbose)

    if out_file: pos_file = "%s.pos" % out_file
    else: pos_file = "vcf.%s.pos" % int(time.time())
    vcf_union_snps(vcfs, pos_file, call_hets, verbose)

    #### Use bcftools to do snarping
    cmd = mpileup_cmd + " -l %s" % pos_file
    for i in range(0, len(names)):
        cmd += " %s" % samples[names[i]]
    cmd += " | bcftools view -"
    if out_file: cmd += " > %s" % out_file
    syscall(cmd, verbose)

    ## remove temporary files
    os.remove(pos_file)
    for filename in vcfs: os.remove(filename)

    return


#####################################################
### Return samtools mpileup cmd
#####################################################
def get_mpileup_cmd(refasta, region, min_qual, vcf, other_params):
    cmd = 'samtools mpileup '
    if vcf: cmd += '-uf %s' % refasta
    else:   cmd += '-f %s' % refasta
    if min_qual: cmd += ' -q %s' % min_qual
    if region:   cmd += ' -r %s' % region
    if other_params: cmd += ' %s' % other_params
    return cmd


##############################
# Main 
##############################
def main():
    options, args = parse_args()
    delimit = '\t'
    refasta = options.refasta
    out_file = options.output
    pos_file = options.pos_file
    vcf = options.vcf
    verbose = options.verbose

    samples = get_samples(options.bam_files, options.bam_list)
    mpileup_cmd = get_mpileup_cmd(refasta, options.region, options.qual, vcf, options.param)

    if pos_file:
        snarp(samples, refasta, pos_file, mpileup_cmd, out_file, delimit, verbose)
    elif vcf:
        variant_call_by_bcftools(samples, mpileup_cmd, out_file, options.call_hets, delimit, verbose)
    else:
        cutoffs = Variant_Cutoffs()
        cutoffs.load_from_options(options)
        variant_call_by_mpileup(samples, mpileup_cmd, out_file, cutoffs, options.call_hets, options.call_indel, delimit, verbose)

    return


if __name__ == '__main__':
    main()

