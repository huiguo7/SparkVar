#!/usr/bin/env python

"""
yzSAM.py
  SAM/BAM parser class.

  Version 1.0 - Read a SAM/BAM file into a SAM class.
  Version 1.1 - moved BasePair class to yzBASE.py
  Version 1.2 - added INDEL calling
  Version 1.3 - added a function to convert mpileup to consensus sequence with indels
  Version 1.4 - added a function to calculate coverage and identity for each hit
                changed to report revcomp if the a hit is reversed
  Version 1.5 - added a function to do snarp

  Last Updated: 04/09/2015
"""

import re, sys, os, string, operator

from yzIUPAC import *
from yzBASE import *
from yzUtil import revcomp


# a class for each alignment line
class Alignment:
    def __init__(self, rec):
        self.read_name = rec[0]
        self.flag = int(rec[1])
        self.refseq = rec[2]
        self.refpos = int(rec[3])
        self.qual = int(rec[4])
        self.CIGAR = rec[5]
        self.mate_refseq = rec[6]
        self.mate_refpos = int(rec[7])
        self.frag_size = int(rec[8])
        self.read_seq = rec[9]
        self.read_qual = rec[10]
        self.optional_fields = string.join(rec[11:],'\t')
        self.options = {}
        for i in range(11, len(rec)):
            val = rec[i].split(':')
            self.options[val[0]] = val[2]

    def edit_dist(self):
        if 'NM' in self.options:
            return int(self.options['NM'])
        return -1

    def mm_code(self):
        if 'MD' in self.options:
            return self.options['MD']
        return None

    def mate_id(self):
        if operator.and_(self.flag, 0x40): return 1
        elif operator.and_(self.flag, 0x80): return 2
        return 0

    def reverse(self):
        if operator.and_(self.flag, 0x10): return True
        return False

    def mm_pos(self):
        if self.edit_dist() == 0: return '-'
        md = self.mm_code()
        if not md: return None
        val = re.split('[ACGTN]+', md, flags=re.IGNORECASE)
        pos = []
        p = 0
        for i in range(0, len(val)-1):
            p += int(val[i]) + 1
            pos.append(p)
        return pos

    def mm_pos_toprint(self):
        rec = (str(i) for i in self.mm_pos())
        return string.join(rec,';')

    def seq_mm_lowercase(self):
        if self.edit_dist() == 0: 
            seq = self.read_seq
        else:
            bases = list(self.read_seq.upper())
            pos = self.mm_pos()
            if not pos:
                seq = self.read_seq
            else:
                for i in xrange(len(pos)):
                    bases[pos[i]-1] = bases[pos[i]-1].lower()
                seq = string.join(bases,'')
        if self.reverse(): seq = revcomp(seq)
        return seq

    def get_read_seq(self):
        seq = self.read_seq
        if self.reverse(): return revcomp(seq)
        else: return seq
       
    def get_read_qual(self):
        qual = self.read_qual
        if self.reverse(): return qual[::-1]
        else: return qual
 
    def aligned(self):
        if not 'AS' in self.options: return False
        return True

    def align_score(self):
        return self.get_read_qual()

    def has_mate(self):
        if operator.mod(self.flag,2) == 1: return True
        return False

    def unique(self):
        if not self.aligned(): return False
        if 'XS' in self.options and 'AS' in self.options:
            if self.options['XS'] == self.options['AS']:
                return False
        return True

    def cigar(self):
        m = re.findall('[0-9]+[MIDNSHP=X]', self.CIGAR)
        return m 

    def set_cov_ident(self):
        rec = self.cigar()
        (total_len, gap_len, match_len) = (0, 0, 0)
        for i in range(0, len(rec)):
            val = rec[i]
            size = int(val[:-1])
            code = val[-1]
            total_len += size
            if code in ['S','H']: gap_len += size
            elif code == 'M': match_len += size
        (self.cov, self.ident) = (0.0, 0.0)
        if match_len > 0: self.ident = 1 - float(self.edit_dist())/float(match_len)
        if total_len > 0: self.cov = 1 - float(gap_len)/float(total_len)
        self.start_pos = self.refpos
        self.stop_pos = self.start_pos + match_len
        return
       
    def duplicate(self, dest, min_dist):
        if self.refseq != dest.refseq: return False
        elif self.mate_refseq != dest.mate_refseq: return False
        else:
            if self.reverse():
                my_pos = self.refpos + len(self.read_seq)
                dest_pos = dest.refpos + len(dest.read_seq)
            else:
                my_pos = self.refpos
                dest_pos = dest.refpos
            if my_pos != dest_pos: return False
            elif self.mate_refpos - dest.mate_refpos > min_dist: return False
            else: return True

 
    def write(self, outfile, delimit='\t'):
        outfile.write("%s" % self.read_name)
        outfile.write("%s%s" % (delimit,self.flag))
        outfile.write("%s%s" % (delimit,self.refseq))
        outfile.write("%s%s" % (delimit,self.refpos))
        outfile.write("%s%s" % (delimit,self.qual))
        outfile.write("%s%s" % (delimit,self.CIGAR))
        outfile.write("%s%s" % (delimit,self.mate_refseq))
        outfile.write("%s%s" % (delimit,self.mate_refpos))
        outfile.write("%s%s" % (delimit,self.frag_size))
        outfile.write("%s%s" % (delimit,self.read_seq))
        outfile.write("%s%s" % (delimit,self.read_qual))
        outfile.write("%s%s" % (delimit,string.join(self.optional_fields.split('\t'),delimit)))
        outfile.write("\n")
        

    def write_info(self, outfile, delimit='\t'):
        outfile.write("%s" % (self.read_name))
        outfile.write("%s%s" % (delimit,self.read_seq))
        outfile.write("%s%s" % (delimit,self.refseq))
        outfile.write("%s%s" % (delimit,self.refpos))
        outfile.write("%s%s" % (delimit,self.edit_dist()))
#        outfile.write("%s%s" % (delimit,self.mm_pos_toprint()))
        outfile.write("\n")


# a class for pair mates
class Pair:
    def __init__(self, a1, a2):
        self.a1 = a1
        self.a2 = a2

    def edit_dist(self):
        return self.a1.edit_dist()+self.a2.edit_dist()

    def write(self, outfile):
        self.a1.write(outfile)
        self.a2.write(outfile)

    def write_info(self, outfile, delimit='\t'):
        outfile.write("%s%s" % (self.a1.read_name,delimit))
        read1_seq = self.a1.seq_mm_lowercase()
        read2_seq = self.a2.seq_mm_lowercase()
        outfile.write("%s%s" % (read1_seq,delimit))
        outfile.write("%s%s" % (read2_seq,delimit))
        outfile.write("%s%s" % (self.a1.refseq,delimit))
        outfile.write("%s%s" % (self.a1.refpos,delimit))
        outfile.write("%s%s" % (self.a1.mate_refpos,delimit))
        outfile.write("%s%s" % (self.a1.frag_size,delimit))
        outfile.write("%s%s" % (self.edit_dist(),delimit))
        outfile.write("%s%s" % (self.a1.edit_dist(),delimit))
        outfile.write("%s" % (self.a2.edit_dist()))
#        outfile.write("%s%s" % (self.a1.mm_pos_toprint(),delimit))
#        outfile.write("%s" % (self.a2.mm_pos_toprint()))
        outfile.write("\n")

    def write_loci(self, outfile, delimit='\t'):
        outfile.write("%s%s" % (self.a1.refseq,delimit))
        outfile.write("%s%s" % (self.a1.refpos,delimit))
        outfile.write("%s%s" % (self.a1.mate_refpos,delimit))
        outfile.write("%s%s" % (self.a1.frag_size,delimit))
        outfile.write("%s%s" % (self.edit_dist(),delimit))
        outfile.write("%s%s" % (self.a1.edit_dist(),delimit))
        outfile.write("%s" % (self.a2.edit_dist()))
        outfile.write("\n")


# a class for each read pairs
class Read:
    def __init__(self, name):
        self.name = name
        self.paired = 0   # by default reads are singles
        self.edit_dist = -1
        self.align = []

    def add_align(self, align):
        self.align.append(align)

    def num_align(self):
        return len(self.align)

    def num_align_mm(self):
        num_pm, num_mm = (0, 0)
        for i in range(0, len(self.align)):
            a = self.align[i]
            if a.edit_dist() == 0: num_pm += 1
            else: num_mm += 1
        return num_pm, num_mm

    def set_as_pair(self):
        self.paired = 1

    def is_pair(self):
        return self.paired

    def best_edit_dist(self):
        self.edit_dist = self.align[0].edit_dist()
        for i in range(1, len(self.align)):
            a = self.align[i]
            if self.edit_dist > a.edit_dist():
                self.edit_dist = a.edit_dist()

    def num_best_align(self):
        if self.edit_dist == -1: self.best_edit_dist()
        num_best = 0
        for i in range(0, len(self.align)):
            a = self.align[i]
            if a.edit_dist() == self.edit_dist:
                num_best += 1
        return num_best

    def best_aligns(self):
        if self.edit_dist == -1: self.best_edit_dist()
        best = []
        for i in range(0, len(self.align)):
            a = self.align[i]
            if a.edit_dist() == self.edit_dist:
                best.append(a)
        return best

    def write_best(self, outfile, num_location):
        best_align = self.best_aligns()
        if self.edit_dist >= 0:
            if num_location < 0 or (num_location > 0 and len(best_align) <= num_location):
                for i in range(0, len(best_align)):
                    best_align[i].write(outfile)
        del best_align

    def write_all(self, outfile):
        for i in xrange(len(self.align)):
            a = self.align[i]
            a.write(outfile)
    
    def write_info_all(self, outfile, delimit='\t'):
        for i in xrange(len(self.align)):
            a = self.align[i]
            a.write_info(outfile,delimit)

    def write_loci(self, outfile, delimit='\t'):
        self.align[0].write_loci(outfile,delimit)
 

# a class for sam file
class Samfile:
    def __init__(self):
        self.headers = []
        self.reads = {}

    def load(self, infile):
        line = infile.readline()
        while line != "":
            if line[0] == '@':
                self.headers.append(line)
            else:
                rec = line[:-1].split('\t')
                read_name = rec[0]
                if not read_name in self.reads:
                    self.reads[read_name] = Read(read_name)
                if rec[5] != '*':  # not empty CIGAR 
                    align = Alignment(rec)
                    if align.has_mate():
                        line = infile.readline()
                        rec2 = line[:-1].split('\t')
                        if rec2[0] != read_name:
                            sys.stderr.write("Warning: Missing mate alignment: %s\n" % read_name)
                            sys.exit(1)
                        align2 = Alignment(rec2)
                        pair = Pair(align, align2)
                        self.reads[read_name].add_align(pair)
                        self.reads[read_name].set_as_pair()
                    else:
                        self.reads[read_name].add_align(align)
            line = infile.readline()
        return


    # write the best matches only
    def write_best(self, outfile, num_location):
        for i in range(0, len(self.headers)):
            outfile.write("%s" % self.headers[i])
        names = self.reads.keys()
        names.sort()
        for read_name in names:
            self.reads[read_name].write_best(outfile,num_location)
        return


    # filter to best matches only and up to the given number of locations
    def filter(self, infile, outfile, num_location):
        prev_read_name = ""
        line = infile.readline()
        while line != "":
            if line[0] == '@':
                outfile.write("%s" % line)
            else:
                rec = line[:-1].split('\t')
                read_name = rec[0]
                align = Alignment(rec)
                if read_name != prev_read_name:
                    if prev_read_name != "": 
                        read.write_best(outfile,num_location)
                    read = Read(read_name)
                    prev_read_name = read_name
                if align.has_mate():
                    line = infile.readline()
                    rec2 = line[:-1].split('\t')
                    if rec2[0] != read_name:
                        sys.stderr.write("Warning: Missing mate alignment: %s\n" % read_name)
                        sys.exit(1)
                    align2 = Alignment(rec2)
                    pair = Pair(align, align2)
                    read.add_align(pair)
                    read.paired = 1
                else:
                    read.add_align(align)
            line = infile.readline()
        read.write_best(outfile,num_location)
        return

    # filter to unique matches only for bowtie2 -M result
    def unique(self, infile, outfile, min_align_score):
        line = infile.readline()
        while line != "":
            if line[0] == '@':
                outfile.write("%s" % line)
            else:
                rec = line[:-1].split('\t')
                align = Alignment(rec)
                if align.unique() and align.align_score() >= min_align_score: 
                    align.write(outfile)
            line = infile.readline()
        return



def SAMpileupLineParser(line):
    rec = line[:-1].split('\t')
    basepair = BasePair()
    basepair.load_ref(rec[0],int(rec[1]),rec[2])
    basepair.load_samples_mpileup(rec[3:])
    return basepair




def getBAMSeq4Region(bam_file, refasta, seqid, region="", min_qual=30, major=True):
    if region == "":
        start = 1
        stop = -1
        r = seqid
    else:
        val = region.split('-')
        start = int(val[0])
        stop  = int(val[1])
        r = "%s:%s" % (seqid, region)
    cmd = "samtools mpileup -q %s -f '%s' -r '%s' %s" % (min_qual, refasta, r, bam_file)
    sys.stderr.write("[yzSAM] getBAMSeq4Region(): %s\n" % cmd)
    infile = os.popen("%s" % cmd)
    seq = "" 
    base_id = start
    line = infile.readline()
    while line != "":
        if line[0] != '#':
            base = SAMpileupLineParser(line)
            if base_id < base.pos:
                for j in range(base_id, base.pos): seq += '.'
                base_id = base.pos
            if major:
                seq += base.major_allele()
            else:
                if base.has_insertion(): seq += '+'
                seq += base.iupac()
            base_id += 1
        line = infile.readline()
    if stop > 0 and base_id < stop:
        for j in range(base_id, stop): seq += '.'
    return seq


def mpileup2consensus(bam_file, refasta, region):
    rec = region.split(':')
    if len(rec) == 1:
        seqid = region
        start = 1
        stop  = -1
    else:
        seqid = string.join(rec[:-1],':')
        val = rec[-1].split('-')
        start = int(val[0])
        stop  = int(val[1])
    min_qual = 30
    cmd = "samtools mpileup -q %s -f '%s' -r '%s' %s" % (min_qual, refasta, region, bam_file)
    sys.stderr.write("[yzSAM] mpileup2consensus(): %s\n" % cmd)
    infile = os.popen("%s" % cmd)
    seq = ""
    base_id = start
    line = infile.readline()
    while line != "":
        if line[0] != '#':
            base = SAMpileupLineParser(line)
            if base_id < base.pos:
                for j in range(base_id, base.pos): seq += '-'
                base_id = base.pos
            seq += base.major_allele()
            seq += base.insertion()
            base_id += 1
        line = infile.readline()
    if stop > 0 and base_id < stop:
        for j in range(base_id, stop): seq += '-'
    return seq


def getBAMSeq4All(bam_file, refasta, seqids, qual, major=True):
    sample_seq = {}
    for id in seqids:
        sample_seq[id] = getBAMSeq4Region(bam_file,refasta,id,'',qual,major)
    return sample_seq


def getBAMCov4Region(bam_file, seqid, region="", min_qual=0):
    if region == "": r = seqid
    else: r += "%s:%s" % (seqid, region)
    cmd = "samtools mpileup -q %s -r '%s' %s" % (min_qual, r, bam_file)
    sys.stderr.write("[yzSAM] getBAMCov4Region(): %s\n" % cmd)
    out_file = os.popen(cmd)
    coverage = {}
    line = out_file.readline()
    while line != "":
        base = SAMpileupLineParser(line)
        if base.cov() > 0:
            coverage[base.pos] = base.cov()
        line = out_file.readline()
    return coverage



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


def SNARP_GetRefbases(pos_file, refasta, verbose):
    refbases = {}
    cmd = 'fasta_extract %s -p -i %s -l 0 | ~zhangy3/tools/parser/fasta/fastaParser.py -a convert' % (refasta, pos_file)
    if verbose: sys.stderr.write('[yzSAM] SNARP_GetRefbases(): %s\n' % cmd)
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


def SNARP_PrintEmptyLine(outfile, seqid, coord, num_samples, refbases, delimit):
    values = [seqid, str(coord), 'SNP']
    pos = '%s:%s' % (seqid,coord)
    if pos in refbases: refbase = refbases[pos]
    else: refbase = 'N'
    values += [refbase]
    for i in xrange(num_samples):
        values += ['.','0','0.0']
    out_line = delimit.join(values)
    outfile.write(out_line+'\n')


def SNARP_HeaderWrite(outfile, names, delimit):
    values = ['SeqID','Position','Type','Refbase']
    for i in range(0, len(names)):
        name = names[i]
        values += ['%s_Basecall'%name,'%s_Cov'%name,'%s_Purity'%name]
    out_line = delimit.join(values)
    outfile.write('#'+out_line+'\n')


def SNARP_BAMs(bamfiles, pos_list, refasta, region, qual, het_cov, het_purity, verbose):
    pos_fn = 'yzSAM-snarp-%s.pos' % os.getpid()
    pos_file = open(pos_fn, 'wb')
    for pos in pos_list:
        rec = pos.split(':')
        pos_file.write('\t'.join(rec)+'\n')
    pos_file.close()
    refbases  = SNARP_GetRefbases(pos_fn, refasta, verbose)
    snarp = {}
    for id in pos_list:
        snarp[id] = {'refbase':refbases[id], 'calls':[]}
    cmd = "samtools mpileup -f %s" % refasta
    if region: cmd += " -r %s" % region
    if qual:   cmd += " -q %s" % qual
    cmd += " -l %s" % pos_fn
    for i in xrange(len(bamfiles)):
        cmd += " %s" % bamfiles[i]
    if verbose: sys.stderr.write('[yzSAM] SNARP_BAMs(): %s\n' % cmd)
    fpt = os.popen(cmd)
    line = fpt.readline()
    while line != '':
        if line[0] != '#':
            basepair = SAMpileupLineParser(line)
            id = '%s:%s' % (basepair.seqid, basepair.pos)
            for i in range(0, len(basepair.samples)):
                sample = basepair.samples[i]
                snarp[id]['calls'].append((sample.basecall(het_cov,het_purity),sample.cov,sample.purity()))
        line = fpt.readline()
    os.remove(pos_fn)
    return snarp

