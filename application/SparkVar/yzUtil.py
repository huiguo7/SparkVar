#!/opt/common/compbio/bin/python

"""
yzUtil.py
  A module for utilities.

"""


import os
import sys
import re
import string
import tempfile
import subprocess
from datetime import datetime

rev_comp_table = string.maketrans("ACGTacgt","TGCAtgca")
def revcomp(seq):
    return string.translate(seq, rev_comp_table)[::-1]


def getCommandOutput(command):
    outfile = tempfile.mktemp( )
    errfile = tempfile.mktemp( )
    cmd = "( %s ) > %s 2> %s" % (command, outfile, errfile)
    sys.stderr.write("%s\n" % cmd)
    os.system(cmd)
    test = file(errfile).readline()
    try:
        if (len(test) != 0):
            raise RuntimeError, '%s' % (file(errfile).read( ))
        return file(outfile).readlines()
    finally:
        os.remove(outfile)
        os.remove(errfile)


def log10(val):
    if val == 0: return -1
    return math.log10(val)


def system_run(cmd, verbose=True):
    if verbose: sys.stderr.write("[%s] #cmd: %s\n" % (str(datetime.now()), cmd))
    os.system(cmd)
    return

def make_dir(dir):
    if not os.path.exists("%s" % dir):
        system_run("mkdir %s" % dir)
    else:
        sys.stderr.write("Warning: folder %s exists.\n" % dir)
    return

def cat_cmd(input_file):
    if re.search(".bz2", input_file):
        cmd = "bzcat %s" % input_file
    elif re.search(".gz", input_file):
        cmd = "gzip -d -c %s" % input_file
    else: cmd = "cat %s" % input_file
    return cmd

def open_file(input_file):
#    if re.search(".bz2", input_file):
    if input_file.endswith('.bz2'):
        cmd = "bzcat %s" % input_file
#    elif re.search(".gz", input_file)
    elif input_file.endswith('.gz'):
        cmd = "gzip -d -c %s" % input_file
    else: cmd = "cat %s" % input_file
    sys.stderr.write("%s\n" % cmd)
    return os.popen(cmd, 'rb')

def open_file_2(input_file):
    args = []
#    if re.search(".bz2", input_file):
    if input_file.endswith('.bz2'):
        args.append('bzip2')
        args.append('-d')
        args.append('-c')
        args.append('-q')
        args.append('-k')
#    elif re.search(".gz", input_file):
    elif input_file.endswith('.gz'):
        args.append('gzip')
        args.append('-d')
        args.append('-c')
    else: args.append('cat')
    args.append(input_file)
    print args
    cat = subprocess.Popen(args,stdout=subprocess.PIPE)
    return cat.stdout

