#!/opt/common/compbio/bin/python

"""
yzBASE.py

  A module for base pair infomration.
  
  Version 1.0 -
  Version 1.3 - added the function for return indels caught in bamfiles
  Version 1.5 - fixed a bug in indel calling
"""

import sys, os, re 
from yzIUPAC import *


ACGT = ['A', 'C', 'G', 'T', 'N', 'a', 'c', 'g', 't', 'n']
MIN_FREQ = 0.1

class Indel:
    def __init__(self, dir):
        self.dir = dir
        self.len = 0
        self.seq = ''

    # return accession of the indel
    def acc(self):
        return self.dir + self.seq


# A class for coverage information for one sample on one base pair
class SampleBase:
    def __init__(self):
        self.cov = 0
        self.major = ''
        self.minor = ''
        self.alleles = {}   # alleles and their coverage, sorted by coverage from high to low
        self.indels = {}

    # load coverage from cdb print line
    def load_cdb(self, rec):
        alt = {}
        for i in range(0, 4):
            cov = int(rec[i])
            self.cov += cov
            if cov > 0: alt[ACGT[i]] = cov
        if self.cov == 0: return
        self.alleles = sorted(alt.iteritems(), key=lambda (k,v):(v,k), reverse=True)

    # load coverage from samtools pileup/mpileup
    def load_pileup(self, rec, refbase):
        self.cov = int(rec[0])
        if self.cov == 0: return
        bases = rec[1]
        quals = rec[2]
        alt = {}
        indels = {}
        i = 0
        while i < len(bases):
            base = bases[i]
            if base == '^':   i += 2
            elif base == '$': i += 1
            elif base == '*': 
                i += 1
                if not '-' in alt: alt['-'] = 1
                else: alt['-'] += 1
            elif base in ['-', '+']:
                indel = Indel(base)
                start = i + 1
                end = start
                while bases[end].isdigit():
                    end += 1
                indel.len = int(bases[start:end])
                indel.seq = bases[end:end+indel.len].upper()
                indel_acc = indel.acc()
                if not indel_acc in indels: indels[indel_acc] = 0
                indels[indel_acc] += 1
                i = end + indel.len
            else:
                if base in ['.', ',']: base = refbase.upper()
                elif base in ACGT:     base = base.upper()
                if base in ACGT:
                    if not base in alt: alt[base] = 1
                    else: alt[base] += 1
                i += 1
        if len(alt) == 0: self.cov = 0
        else: self.alleles = sorted(alt.iteritems(), key=lambda (k,v):(v,k), reverse=True)
        if self.cov > 0:
            self.major = self.alleles[0][0].upper()
            if self.num_alleles() > 1:
                if self.alleles[0][1] == self.alleles[1][1] and self.major == refbase.upper():
                    self.minor = self.major
                    self.major = self.alleles[1][0].upper()
                else:
                    self.minor = self.alleles[1][0].upper()
        if len(indels) > 0: 
            self.indels = sorted(indels.iteritems(), key=lambda (k,v):(v,k), reverse=True)

    def covered(self):
        if self.cov == 0: return False
        else: return True

    def purity(self):
        if self.cov == 0: return 0
        return float(self.alleles[0][1])/float(self.cov)

    def ref_freq(self, refbase):
        if self.cov == 0: return 0
        for i in range(0, self.num_alleles()):
            if self.alleles[i][0].upper() == refbase.upper():
                return float(self.alleles[i][1])/float(self.cov)
        return 0

    def alt_freq(self, refbase):
        if self.cov == 0: return 0
        for i in range(0, self.num_alleles()):
            if self.alleles[i][0].upper() != refbase.upper():
                return float(self.alleles[i][1])/float(self.cov)
        return 0

    def num_alleles(self):
        return len(self.alleles)
 
    def num_indels(self):
        return len(self.indels)

    def has_indel(self):
        if len(self.indels) > 0: return True
        return False

    def indel_purity(self):
        if not self.has_indel(): return 0
        indel_cov = 0
        for i in range(0, len(self.indels)):
            indel_cov += self.indels[i][1]
        return float(indel_cov)/float(self.cov)

    def major_allele(self):
        if self.cov == 0: return '.'
        return self.major

    def minor_allele(self):
        return self.minor

    def has_insertion(self):
        if not self.has_indel(): return False
        if self.indels[0][0][0] == '+': return True
        return False

    def insertion(self):
        if self.has_insertion(): return self.indels[0][0][1:]
        return ''
        
    def basecall(self, het_cov, het_purity):
        if self.cov == 0: return '.'
        if self.cov >= het_cov and self.purity() <= het_purity:
#            minor_allele = self.alleles[1][0]
#            minor_allele_perc = self.alleles[1][1]
            return "%s/%s" % (self.major, self.minor)
        else:
            return self.major

    def basecall_all(self):
        if self.cov == 0: return '.'
        basecall = self.alleles[0][0]
        for i in range(1, len(self.alleles)):
            if float(self.alleles[i][1])/float(self.cov) >= MIN_FREQ:
                basecall += "/%s" % self.alleles[i][0]
        return basecall

    def iupac(self):
        basecall = self.basecall_all()
        if len(basecall) == 1: return basecall
        alleles = basecall.split('/')
        return IUPACcode(alleles[0],alleles[1])

#        if self.cov == 0: return '.'
#        if self.cov >= het_cov and self.purity() <= het_purity: 
#            if self.num_alleles() == 1:
#                return self.alleles[0][0]
#            return IUPACcode(self.alleles[0][0],self.alleles[1][0])
#        return self.alleles[0][0]

    def is_snp(self, refbase, cov, purity):
        if self.cov == 0: return False
        if refbase.upper() == 'N': return False
        major_allele = self.major
        if major_allele == 'N' or major_allele == '-': return False
        if self.num_alleles() > 1 and self.alleles[1][0].upper() == '-': return False
        if major_allele != refbase.upper() and self.cov >= cov and self.purity() >= purity:
            return True
        return False

    def is_hets(self, cov, purity):
        if self.cov == 0: return False
        major_allele = self.major
        if major_allele == 'N' or major_allele == '-': return False
        if self.num_alleles() > 1 and self.alleles[1][0].upper() == '-': return False
        if self.cov >= cov and self.purity() <= purity: return True
        return False
    
    def is_varient(self, refbase, cov, max_ref_freq):
        if self.cov == 0: return False
        if self.alleles[0][0].upper() == 'N': return False
        if self.cov >= cov and self.ref_freq(refbase) <= max_ref_freq: return True
        return False
    
    def is_indel(self, cov, purity):
        if self.cov == 0: return False
        if len(self.indels) == 0: return False
        if self.cov >= cov and self.indel_purity() >= purity: return True
        return False

    def max_indel(self):
        if not self.has_indel(): return ''
        max_indel = ''
        num_del = 0
        for i in range(0, len(self.indels)):
            indel = self.indels[i][0]
            if indel[0] == '-' and len(indel)-1 > len(max_indel):
                max_del = indel[1:]
                num_del += 1    
        if len(self.indels) > 0 and num_del/len(self.indels) > 0.5:
            max_indel = max_del
        return max_indel
            
    def write(self, outfile, het_cov, het_purity):
        outfile.write("\t%s" % self.basecall(het_cov,het_purity))
        outfile.write("\t%s" % self.cov)
        outfile.write("\t%.4f" % self.purity())

    def write_indel(self, outfile, refbase, max_indel, indel_perc):
        if self.has_indel(): 
            indel = self.indels[0][0]
            outfile.write("\t%s" % refbase)
            if indel[0] == '+': outfile.write("%s" % indel[1:])
            else: 
                start = len(indel[1:])
                outfile.write("%s" % max_indel[start:])
            for i in range(1, len(self.indels)):
                if (self.indels[i][1]/self.cov) >= indel_perc:
                    indel = self.indels[i][0]
                    outfile.write(",%s" % refbase)
                    if indel[0] == '+': outfile.write("%s" % indel[1:])
                    else: 
                        start = len(indel[1:])
                        outfile.write("%s" % max_indel[start:])
            outfile.write("\t%s" % self.cov)
            outfile.write("\t%.4f" % self.indel_purity()) 
        else:
            outfile.write("\t%s%s" % (refbase,max_indel))
            outfile.write("\t%s" % self.cov)
            outfile.write("\t%.4f" % self.purity())

    def write_varient(self, outfile, refbase):
        outfile.write("\t%s" % self.basecall_all())
        outfile.write("\t%s" % self.cov)
        outfile.write("\t%.4f" % self.ref_freq(refbase))
        acgt = {'A':0, 'C':0, 'G':0, 'T':0}
        for i in range(0, self.num_alleles()):
            acgt[self.alleles[i][0]] = self.alleles[i][1]
        outfile.write("\t%s,%s,%s,%s" % (acgt['A'],acgt['C'],acgt['G'],acgt['T']))


# A class for one base pair, including reference information and
# coverage infomration for a set of samples
class BasePair:
    def __init__(self):
        self.seqid = ""
        self.pos = -1
        self.refbase = '.'
        self.samples = []

    def load_ref(self, id, pos, ref):
        self.seqid = id
        self.pos = pos
        self.refbase = ref

    def load_samples_mpileup(self, rec):
        for i in range(0, len(rec), 3):
            base = SampleBase()
            base.load_pileup(rec[i:i+3],self.refbase)
            self.samples.append(base)

    # input is the four-tuple of coverage
    def load_samples_cdb(self, rec):
        base = SampleBase()
        base.load_cdb(rec)
        self.samples.append(base)

    def num_samples(self):
        return len(self.samples)

    def is_snp(self, cov, purity):
        if self.refbase == '.' or self.refbase == 'N': return False
        for i in range(0, len(self.samples)):
            if self.samples[i].is_snp(self.refbase,cov,purity): return True
        return False

    def is_hets(self, cov, purity):
        if self.refbase == '.' or self.refbase == 'N': return False
        for i in range(0, len(self.samples)):
            if self.samples[i].is_hets(cov,purity): return True
        return False

    def is_indel(self, cov, purity):
        if self.refbase == '.' or self.refbase == 'N': return False
        for i in range(0, len(self.samples)):
            if self.samples[i].is_indel(cov, purity): return True
        return False

    def is_varient(self, cov, purity):
        if self.refbase == '.' or self.refbase == 'N': return False
        for i in range(0, len(self.samples)):
            if self.samples[i].is_varient(self.refbase,cov,purity): return True
        return False
        
    def covered(self, sample_id = 0):
        if self.num_samples() == 0: return False
        return self.samples[sample_id].covered()

    def cov(self, sample_id = 0):
        if len(self.samples) == 0: return 0
        return self.samples[sample_id].cov

    def basecall(self, sample_id, het_cov, het_purity):
        if self.num_samples() == 0: return '.'
        return self.samples[sample_id].basecall(het_cov,het_purity)

    def purity(self, sample_id = 0):
        if self.num_samples() == 0: return -1
        return self.samples[sample_id].purity()

    def iupac(self, sample_id = 0):
        if self.num_samples() == 0: return '.'
        return self.samples[sample_id].iupac()

    def major_allele(self, sample_id = 0):
        if self.num_samples() == 0: return '.'
        return self.samples[sample_id].major_allele()

    def minor_allele(self, sample_id = 0):
        if self.num_samples() == 0: return ''
        return self.samples[sample_id].minor_allele()

    def has_insertion(self, sample_id = 0):
        if self.num_samples() == 0: return False
        return self.samples[sample_id].has_insertion()

    def insertion(self, sample_id = 0):
        if self.num_samples() == 0: return ''
        return self.samples[sample_id].insertion()

    def indel_purity(self, sample_id = 0):
        if self.num_samples() == 0: return ''
        return self.samples[sample_id].indel_purity()

    def max_indel(self):
        max_indel = ''
        (num_del, num_ins) = (0, 0)
        for i in range(0, len(self.samples)):
            if self.samples[i].has_indel():
                max = self.samples[i].max_indel()
                if len(max) > 0: 
                    num_del += 1
                    if len(max) > len(max_indel): max_indel = max
                else:
                    num_ins += 1
        if num_del > 0 and num_ins > 0: type = 'IND'
        elif num_del > 0: type = 'DEL'
        elif num_ins > 0: type = 'INS'
        return (type, max_indel)

    def write_pos(self, outfile):
        outfile.write("%s\t%s\t" % (self.seqid, self.pos))

    def write(self, outfile, het_cov, het_purity):
        self.write_snp(outfile, het_cov, het_purity)

    def write_snp(self, outfile, het_cov, het_purity):
        self.write_pos(outfile)
        outfile.write("SNP\t")
        outfile.write("%s" % self.refbase)
        for i in range(0, len(self.samples)):
            self.samples[i].write(outfile, het_cov, het_purity)
        outfile.write("\n")

    def write_hets(self, outfile, het_cov, het_purity):
        self.write_pos(outfile)
        outfile.write("HET\t")
        outfile.write("%s" % self.refbase)
        for i in range(0, len(self.samples)):
            self.samples[i].write(outfile, het_cov, het_purity)
        outfile.write("\n")

    def write_indel(self, outfile, indel_perc):
        (type, max_indel) = self.max_indel()
        self.write_pos(outfile)
        outfile.write('%s\t' % type)
        outfile.write("%s%s" % (self.refbase,max_indel))
        for i in range(0, len(self.samples)):
            self.samples[i].write_indel(outfile,self.refbase,max_indel,indel_perc)
        outfile.write("\n")
 
    def write_varient(self, outfile):
        self.write_pos(outfile)
        outfile.write("BSA\t")
        outfile.write("%s" % self.refbase)
        for i in range(0, len(self.samples)):
            self.samples[i].write_varient(outfile,self.refbase)
        outfile.write("\n")

