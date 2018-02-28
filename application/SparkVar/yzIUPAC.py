#!/opt/common/compbio/bin/python

import os, sys, re, string

# IUPAC ambiguous nucleotide code
#A       A       Adenine
#C       C       Cytosine
#G       G       Guanine
#T       T       Thymine
#R       A or G  puRine
#Y       C or T  pYrimidine
#W       A or T  Weak hydrogen bonding
#S       C or G  Strong hydrogen bonding
#M       A or C  aMino group at common position
#K       G or T  Keto group at common position
#B       C, G or T       not A
#D       A, G or T       not C
#H       A, C or T       not G
#V       A, C or G       not T
#N       A, C, G or T    aNy

DeIUPAC = {'A':'A','C':'C','G':'G','T':'T',
           'R':['A','G'],'Y':['C','T'],'W':['A','T'],
           'S':['C','G'],'M':['A','C'],'K':['G','T'],
           'B':['C','G','T'],'D':['A','G','T'],
           'H':['A','C','T'],'V':['A','C','G'],
           'N':['A','C','G','T']}

UniIUPAC = {'[A/A]':'A','[C/C]':'C','[G/G]':'G','[T/T]':'T',
            '[A/-]':'A','[C/-]':'C','[G/-]':'G','[T/-]':'T',
            '[A/:]':'A','[C/:]':'C','[G/:]':'G','[T/:]':'T',
            '[:/A]':'A','[:/C]':'C','[:/G]':'G','[:/T]':'T',
            '[-/A]':'A','[-/C]':'C','[-/G]':'G','[-/T]':'T',
            '[A]':'A','[C]':'C','[G]':'G','[T]':'T',
            '[A/]':'A','[C/]':'C','[G/]':'G','[T/]':'T'}

BiIUPAC = {'[A/G]':'R','[C/T]':'Y','[A/T]':'W',
           '[G/A]':'R','[T/C]':'Y','[T/A]':'W',
           '[C/G]':'S','[A/C]':'M','[G/T]':'K',
           '[G/C]':'S','[C/A]':'M','[T/G]':'K',
           '[AG]':'R','[CT]':'Y','[AT]':'W',
           '[GA]':'R','[TC]':'Y','[TA]':'W',
           '[CG]':'S','[AC]':'M','[GT]':'K',
           '[GC]':'S','[CA]':'M','[TG]':'K'}

TriIUPAC = {'[C/G/T]':'B','[C/T/G]':'B','[G/C/T]':'B',
         '[G/T/C]':'B','[T/C/G]':'B','[T/G/C]':'B',
         '[A/G/T]':'D','[A/T/G]':'D','[G/A/T]':'D',
         '[G/T/A]':'D','[T/A/G]':'D','[T/G/A]':'D',
         '[A/C/T]':'H','[A/T/C]':'H','[C/A/T]':'H',
         '[C/T/A]':'H','[T/A/C]':'H','[T/C/A]':'H',
         '[A/G/C]':'V','[A/C/G]':'V','[G/A/C]':'V',
         '[G/C/A]':'V','[C/A/G]':'V','[C/G/A]':'V',
         '[CGT]':'B','[CTG]':'B','[GCT]':'B',
         '[GTC]':'B','[TCG]':'B','[TGC]':'B',
         '[AGT]':'D','[ATG]':'D','[GAT]':'D',
         '[GTA]':'D','[TAG]':'D','[TGA]':'D',
         '[ACT]':'H','[ATC]':'H','[CAT]':'H',
         '[CTA]':'H','[TAC]':'H','[TCA]':'H',
         '[AGC]':'V','[ACG]':'V','[GAC]':'V',
         '[GCA]':'V','[CAG]':'V','[CGA]':'V'}

IUPAC2Bracket = {'A':'A','C':'C','G':'G','T':'T',
         'R':"[A/G]",'Y':"[C/T]",'W':"[A/T]",
         'S':"[C/G]",'M':"[A/C]",'K':"[G/T]",
         'B':"[C/G/T]",'D':"[A/G/T]",
         'H':"[A/C/T]",'V':"[A/C/G]",
         'N':"[A/C/G/T]"}

def SeqConvert2IUPAC(in_seq, verbose=True):
    seq = in_seq.upper()
    # added for cases like [A/N]
    if re.search("\[[ACGT]\/N\]", seq):
        return re.sub(r'\[[ACGT]\/N\]','N',seq)
    m = re.search("\[([ACGT-])*(\/)*([ACGT-])*\]", seq)
    if not m: 
        if verbose: sys.stderr.write("[yzIUPAC] Warning: invalid SNP: %s\n" % seq)
        return seq
    snp_call = m.group(0)
#    alleles = snp_call[1:-1].split('/')
    if snp_call in UniIUPAC:
        iupac = UniIUPAC[snp_call]
    elif snp_call in BiIUPAC:
        iupac = BiIUPAC[snp_call]
    elif snp_call in TriIUPAC:
        iupac = TriIUPAC[snp_call]
    else:
        if verbose: sys.stderr.write("[yzIUPAC] Warning: invalid SNP: %s\n" % snp_call)
        return seq
    new_seq = re.sub(r'\[[ACGT-]((\/)*[ACGT-])*\]', iupac, seq)
    return new_seq

# input sequence should have SNP in bracket
def SNPosition(seq):
    val = seq.split('[')
    if len(val) > 1:
        return len(val[0]) + 1
    else:
        return len(val[0])/2

# check to see if input sequence has IUPAC in it
def SeqHasIUPAC(seq):
    upper_seq = seq.upper()
    for i in range(0, len(upper_seq)):
        base = upper_seq[i]
        if not base in UniIUPAC.values():
            if base in TriIUPAC.values() or base in BiIUPAC.values(): return True
            else:
                sys.stderr.write("[yzIUPAC] Error: Invalid code in sequence %s.\n" % seq)
                sys.exit(1)
    return False


def IUPACcode(al1, al2):
    allele1 = string.upper(al1)
    allele2 = string.upper(al2)
    unknown = ['-','.','U']
    if (allele1 in unknown) and allele2 in DeIUPAC.keys(): return allele2
    if (allele2 in unknown) and allele1 in DeIUPAC.keys(): return allele1
    if (allele1 in unknown) and (allele2 in unknown):
        if allele1 == '-' or allele2 == '-': return '-'
        if allele1 == 'U' or allele2 == 'U': return 'U'
        return '.'
    if not allele1 in DeIUPAC.keys() or not allele2 in DeIUPAC.keys():
        sys.stderr.write("[yzIUPAC] Warning: Unexpected code found in %s %s\n" % (allele1, allele2))
        return 'U'
    if allele1 == 'A' and allele2 == 'A': return 'A'
    if allele1 == 'C' and allele2 == 'C': return 'C'
    if allele1 == 'G' and allele2 == 'G': return 'G'
    if allele1 == 'T' and allele2 == 'T': return 'T'
    if (allele1 == 'A' and allele2 == 'G') or (allele1 == 'G' and allele2 == 'A'): return 'R'
    if (allele1 == 'C' and allele2 == 'T') or (allele1 == 'T' and allele2 == 'C'): return 'Y'
    if (allele1 == 'A' and allele2 == 'T') or (allele1 == 'T' and allele2 == 'A'): return 'W'
    if (allele1 == 'C' and allele2 == 'G') or (allele1 == 'G' and allele2 == 'C'): return 'S'
    if (allele1 == 'A' and allele2 == 'C') or (allele1 == 'C' and allele2 == 'A'): return 'M'
    if (allele1 == 'G' and allele2 == 'T') or (allele1 == 'T' and allele2 == 'G'): return 'K'
    if (allele1 == 'R' and allele2 == 'A') or (allele1 == 'A' and allele2 == 'R') or \
       (allele1 == 'R' and allele2 == 'G') or (allele1 == 'G' and allele2 == 'R'): return 'R'
    if (allele1 == 'Y' and allele2 == 'C') or (allele1 == 'C' and allele2 == 'Y') or \
       (allele1 == 'Y' and allele2 == 'T') or (allele1 == 'T' and allele2 == 'Y'): return 'Y'
    if (allele1 == 'W' and allele2 == 'A') or (allele1 == 'A' and allele2 == 'W') or \
       (allele1 == 'W' and allele2 == 'T') or (allele1 == 'T' and allele2 == 'W'): return 'W'
    if (allele1 == 'S' and allele2 == 'C') or (allele1 == 'C' and allele2 == 'S') or \
       (allele1 == 'S' and allele2 == 'G') or (allele1 == 'G' and allele2 == 'S'): return 'S'
    if (allele1 == 'M' and allele2 == 'A') or (allele1 == 'A' and allele2 == 'M') or \
       (allele1 == 'M' and allele2 == 'C') or (allele1 == 'C' and allele2 == 'M'): return 'M'
    if (allele1 == 'K' and allele2 == 'G') or (allele1 == 'G' and allele2 == 'K') or \
       (allele1 == 'K' and allele2 == 'T') or (allele1 == 'T' and allele2 == 'K'): return 'K'
    if (allele1 == 'S' and allele2 == 'T') or (allele1 == 'T' and allele2 == 'S') or \
       (allele1 == 'K' and allele2 == 'C') or (allele1 == 'C' and allele2 == 'K') or \
       (allele1 == 'Y' and allele2 == 'G') or (allele1 == 'G' and allele2 == 'Y'): return 'B'
    if (allele1 == 'R' and allele2 == 'T') or (allele1 == 'T' and allele2 == 'R') or \
       (allele1 == 'K' and allele2 == 'A') or (allele1 == 'A' and allele2 == 'K') or \
       (allele1 == 'W' and allele2 == 'G') or (allele1 == 'G' and allele2 == 'W'): return 'D'
    if (allele1 == 'M' and allele2 == 'T') or (allele1 == 'T' and allele2 == 'M') or \
       (allele1 == 'Y' and allele2 == 'A') or (allele1 == 'A' and allele2 == 'Y') or \
       (allele1 == 'W' and allele2 == 'C') or (allele1 == 'C' and allele2 == 'W'): return 'H'
    if (allele1 == 'M' and allele2 == 'G') or (allele1 == 'G' and allele2 == 'M') or \
       (allele1 == 'S' and allele2 == 'A') or (allele1 == 'A' and allele2 == 'S') or \
       (allele1 == 'R' and allele2 == 'C') or (allele1 == 'C' and allele2 == 'R'): return 'V'
    if (allele1 == 'B' and allele2 == 'T') or (allele1 == 'T' and allele2 == 'B') or \
       (allele1 == 'B' and allele2 == 'C') or (allele1 == 'C' and allele2 == 'B') or \
       (allele1 == 'B' and allele2 == 'G') or (allele1 == 'G' and allele2 == 'B'): return 'B'
    if (allele1 == 'D' and allele2 == 'T') or (allele1 == 'T' and allele2 == 'D') or \
       (allele1 == 'D' and allele2 == 'A') or (allele1 == 'A' and allele2 == 'D') or \
       (allele1 == 'D' and allele2 == 'G') or (allele1 == 'G' and allele2 == 'D'): return 'D'
    if (allele1 == 'H' and allele2 == 'T') or (allele1 == 'T' and allele2 == 'H') or \
       (allele1 == 'H' and allele2 == 'A') or (allele1 == 'A' and allele2 == 'H') or \
       (allele1 == 'H' and allele2 == 'C') or (allele1 == 'C' and allele2 == 'H'): return 'H'
    if (allele1 == 'V' and allele2 == 'G') or (allele1 == 'G' and allele2 == 'V') or \
       (allele1 == 'V' and allele2 == 'A') or (allele1 == 'A' and allele2 == 'V') or \
       (allele1 == 'V' and allele2 == 'C') or (allele1 == 'C' and allele2 == 'V'): return 'V'
    return 'N'


def Set2IUPAC(bases):
    iupac = bases[0]
    for i in range(1, len(bases)):
        iupac = IUPACcode(bases[i], iupac)
    return iupac


# get consensus of two sequence
def ConsensusSeq(seq1, seq2):
    if seq1 == "": return seq2
    if seq2 == "": return seq1
    seq = []
    if len(seq1) > len(seq2): min_len = len(seq2)
    else: min_len = len(seq1)
    for i in range(0, min_len):
        if seq2[i] != seq1[i] and seq2[i] != ".":
            basecall = IUPACcode(seq1[i],seq2[i])
        else:
            basecall = seq1[i]
        seq.append(basecall)
    if len(seq1) > len(seq2):
        for i in range(len(seq2),len(seq1)):
            seq.append(seq1[i])
    return string.join(seq,'')

