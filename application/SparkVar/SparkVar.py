"""
SparkVar utlizes spark programming model to speed up reads alignment and variant calling on AWS EMR cluster.
"""

import json
import math
import multiprocessing
import argparse
import os
from pyspark import SparkConf, SparkContext
import re
import sys
import time
import bisect
import functools
from subprocess import call, check_output, Popen, PIPE

# Get arguments from the command line
def parse_args():
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                     description='''\
Welcome to
      ____              __  __    __   
     / __/__  ___  ____/ /__\ '  / /__  ____
    _\ \/ _ \/ _ `/ __/  '_/ \ '/ / _ '/ __/
   /__ / .__/\_,_/_/ /_/\_\   \__/\_,_/_/
      /_/                                          ''')

    Defaults = {'MinMapQual':30, 'MinCov':3, 'MinPurity':0.98, 'HetsCov':10, 'HetsPurity':0.7, \
                'IndelCov':1, 'IndelPurity':0.5, 'IndelPerc':0.1}

    parser.add_argument("--sample", metavar='STR', required=True,
        help="required input sample s3 path (multiple paths delimited by ',')")
    parser.add_argument("--bt2", metavar='STR', required=True,
        help="required reference bowtie2 index prefix")
    parser.add_argument("--ref", metavar='STR', required=True,
        help="required reference sequence file")
    parser.add_argument("--output", metavar='STR', required=True,
        help="required s3 object path for output")
    parser.add_argument("--job_id", metavar='STR', default='SparkVar_job',
        help="optional job id")
    parser.add_argument("--lib_type", metavar='STR', default='paired_end', choices=['single_end', 'paired_end'],
        help="Type of DNA library, choices are [single_end, paired_end]. [paired_end]")
    parser.add_argument("--bt2_param", metavar='\'STR\'',
        help="option to set bowtie2 parameters.")
    parser.add_argument("--save_BAMs", default=False, action="store_true",
        help="Flag to save BAM files after processing is complete. [False]")
    parser.add_argument("--index", default=False, action="store_true",
        help="Option to index the final bam file. [False]")
    parser.add_argument("--aln_only", default=False, action="store_true",
        help="Option to only perform alignment (force '--save_BAMs'). [False]")
    parser.add_argument("--number_partitions", metavar='INT', type=int,
        help="Number of data partitions. [inferred from type of instances and amount of data]")
    parser.add_argument("--qual", metavar='INT', type=int, default=Defaults['MinMapQual'],
        help="Option to set the minimum mapping quality for SNP call. [{}]".format(Defaults['MinMapQual']))
    parser.add_argument("--indel", default=False, action='store_true',
        help="Flag to call INDELs too, defaults to not call indels. [False]")
    parser.add_argument("--hets", default=False, action='store_true',
        help="Flag to call heterozygous SNPs. [False]")
    parser.add_argument("--vcf", default=False, action='store_true',
        help="Option to output vcf files. [tabular]")
    parser.add_argument("--mpileup_param", metavar='\'STR\'',
        help="Option to set samtools mpileup parameters.")
    parser.add_argument("--min_cov", metavar='INT', default=Defaults['MinCov'], type=int,
        help="Option to set the minimum coverage for snp call. [{}]".format(Defaults['MinCov']))
    parser.add_argument("--min_purity", metavar='FLOAT', default=Defaults['MinPurity'], type=float,
        help="Option to set the minimum purity for snp call. [{}]".format(Defaults['MinPurity']))
    parser.add_argument("--hets_cov", metavar='INT', default=Defaults['HetsCov'], type=int,
        help="Option to set the minimum coverage for hets call. [{}]".format(Defaults['HetsCov']))
    parser.add_argument("--hets_purity", metavar='FLOAT', default=Defaults['HetsPurity'], type=float,
        help="Option to set the maximum purity for hets call. [{}]".format(Defaults['HetsPurity']))
    parser.add_argument("--indel_cov", metavar='INT', default=Defaults['IndelCov'], type=int,
        help="Option to set the minimum coverage for indel call. [{}]".format(Defaults['IndelCov']))
    parser.add_argument("--indel_purity", metavar='FLOAT', default=Defaults['IndelPurity'], type=float,
        help="Option to set the minimum purity for indel call. [{}]".format(Defaults['IndelPurity']))
    parser.add_argument("--indel_perc", metavar='FLOAT', default=Defaults['IndelPerc'], type=float,
        help="Option to set the minimum percent for displaying an indel allele. [{}]".format(Defaults['IndelPerc']))
   
    args = vars(parser.parse_args())

    return args

# Convert the AVRO file to an RDD
def get_flatq_rdd(spark_context, path):
    """
    SparkContext -> String -> RDD FLATQ

    Retreives an RDD of FLATQ records. 'path' is the path to input file(s).
    Accepts comma separated or wildcards to specify multiple input paths.
    """
    return spark_context.newAPIHadoopFile(path=path,
            inputFormatClass="org.apache.avro.mapreduce.AvroKeyInputFormat",
            keyClass="org.apache.avro.mapred.AvroKey",
            valueClass="org.apache.hadoop.io.NullWritable",
            keyConverter="org.apache.spark.examples.pythonconverters.AvroWrapperToJavaConverter")

def get_pair_rdd(sc, paths):
    """
    Create pair rdd (filename, record)
    """
    def proc(f):
        gt = os.path.basename(f)[:-5]
        return get_flatq_rdd(sc, f).map(lambda x: (gt, x))

    rdd = functools.reduce(
        lambda rdd1, rdd2: rdd1.union(rdd2),
        (proc(f) for f in paths))
    return rdd

def get_threads(executor_memory, executor_cores):
    # Get number of vCPUs and system memory
    with open('/mnt/var/lib/info/job-flow.json', 'r') as f:
        config_data = json.load(f)
        instance_type = config_data['slaveInstanceType']

    # Memory availbe to Spark for different EC2 instance types
    instance_types = {
        'm1.xlarge': {'vcpu': 4, 'memory': 12288},
        'cc2.8xlarge': {'vcpu': 32, 'memory': 56320},
        'c3.8xlarge': {'vcpu': 32, 'memory': 53248},
        'c4.8xlarge': {'vcpu': 36, 'memory': 53248},
        'm4.10xlarge': {'vcpu': 40, 'memory': 155648},
        'm4.16xlarge': {'vcpu': 64, 'memory': 253952},
    }

    core_count = multiprocessing.cpu_count()
    ram_value = float(instance_types[instance_type]['memory']) / 1024
    ram_value -= 0.1 * ram_value

    # Calculate the nubmer of threads available to the worker
    threads = int(math.floor(core_count / (executor_cores * math.floor(ram_value / executor_memory))))
    return threads

def aln_tab5(seqi, seqs, args, executor_memory, executor_cores, gt):
    # get number of threads for bowtie2
    threads = get_threads(executor_memory, executor_cores)
    print('\n\n\tData are being processed as {:,} threads\n\n'.format(threads))

    # bowtie2 index
    ref =  '/mnt/{}'.format(args['bt2'])

    # Build Bowtie2 command
    bt2_cmd = "bowtie2 -p {} -x {} --quiet --no-hd --tab5 - ".format(threads, ref)

    # Add any user specified options
    if args['bt2_param']: bt2_cmd += " {}".format(args['bt2_param'])

    bt2_proc = Popen(bt2_cmd, stdin=PIPE, stdout=PIPE, bufsize=0, shell=True)

    for seq in seqs:
        seq = seq[0]
        read = seq['header']+'\t'+seq["seq1"]+'\t'+seq["qual1"]+'\t'+seq["seq2"]+'\t'+seq["qual2"]+'\n'
        bt2_proc.stdin.write(bytes(read, 'UTF-8'))

    sam_records = bt2_proc.communicate()[0].strip().decode("utf-8")

    bt2_proc.stdin.close()
    bt2_proc.stdout.close()

    yield (gt, sam_records)

def aln(seqi, seqs, args, executor_memory, executor_cores, gt):
    # get number of threads for bowtie2
    threads = get_threads(executor_memory, executor_cores)
    print('\n\n\tData are being processed as {:,} threads\n\n'.format(threads))

    # bowtie2 index
    ref =  '/mnt/{}'.format(args['bt2'])
    
    # Build Bowtie2 command
    if args['lib_type'] == 'paired_end':
        bt2_cmd = "bowtie2 -p {} -x {} --quiet --no-hd --interleaved - ".format(threads, ref)
    else:
        bt2_cmd = "bowtie2 -p {} -x {} --quiet --no-hd -U - ".format(threads, ref)
    
    # Add any user specified options
    if args['bt2_param']: bt2_cmd += " {}".format(args['bt2_param'])
    
    bt2_proc = Popen(bt2_cmd, stdin=PIPE, stdout=PIPE, bufsize=0, shell=True)

    for seq in seqs:
        seq = seq[0]
        if args['lib_type'] == 'paired_end':
            read = '@'+seq['header']+'_1\n'+seq["seq1"]+'\n+\n'+seq["qual1"]+'\n'+ \
                   '@'+seq['header']+'_2\n'+seq["seq2"]+'\n+\n'+seq["qual2"]+'\n'
        else:
            read = '@'+seq['header']+'\n'+seq["seq1"]+'\n+\n'+seq["qual1"]+'\n'
        bt2_proc.stdin.write(bytes(read, 'UTF-8'))

    sam_records = bt2_proc.communicate()[0].strip().decode("utf-8")

    bt2_proc.stdin.close()
    bt2_proc.stdout.close()

    yield (gt, sam_records)

def split_sam_records(x):
    k, v = x
    ret = []
    for i in v.split('\n'):
        ret.append((k,i))
    return ret

def reads2aln(sc, sample, args, min_partitions, executor_memory, executor_cores):
    # get genotype id
    gt = os.path.basename(sample)[:-5]

    # Caluclate number of partitions based on size of input files
    input_size = float(check_output('aws s3 ls {}'.format(sample), shell=True, universal_newlines=True).split()[2])

    # Repartition based on file size, with a ceiling of 2500
    partitions = math.ceil(input_size/60000)
    partitions = max(partitions, min_partitions)
    if args['number_partitions']: partitions = args['number_partitions']
    print('\n\n\tSample is being processed as {:,} partitions\n\n'.format(partitions))

    # get data into RDD
    step_time = time.time()
    fq_rdd = get_flatq_rdd(sc, sample)

    # Repartition the data so that we don't have more than 4000 partitions (Samtools merge limit)
    fq_rdd = fq_rdd.repartition(partitions)
    print("---create rdd from avro {} seconds ---".format(time.time() - step_time))

    # Execute Bowtie2 alignment
    step_time = time.time()
    sam_rdd = fq_rdd.mapPartitionsWithIndex(lambda seqi, seqs: aln_tab5(seqi, seqs, args, executor_memory, executor_cores, gt))\
                    .flatMap(split_sam_records)
    print("---map step {} seconds ---".format(time.time() - step_time))

    return sam_rdd, partitions

def get_sam_header(args):
    ref =  '/mnt/{}'.format(args['bt2'])
    if args['lib_type'] == 'paired_end':
        bt2_cmd = "bowtie2 --quiet -x {} --interleaved - ".format(ref)
        read = '@seq_1\nACGT\n+\nIIII\n@seq_2\nACGT\n+\nIIII\n'
    else:
        bt2_cmd = "bowtie2 -x {} --no-hd -U - ".format(ref)
        read = '@fqid\nACGT\n+\nIIII\n'
    if args['bt2_param']: bt2_cmd += " {}".format(args['bt2_param'])
    bt2_proc = Popen(bt2_cmd, stdin=PIPE, stdout=PIPE, bufsize=0, shell=True)
    bt2_proc.stdin.write(bytes(read, 'UTF-8'))
    sam_records = bt2_proc.communicate()[0].strip().decode("utf-8")
    bt2_proc.stdin.close()
    bt2_proc.stdout.close()
    header = [i for i in sam_records.split('\n') if i.startswith('@')]
    return '\n'.join(header)

def get_chrom_info(header):
    def sort_chroms(x):
        for i in x.split('\t'):
            if i.startswith('SN'):
                return int(re.search('\d+', i[3:]).group())

    header_lines = header.split('\n')
    chrom_lines = [i for i in header_lines if i.startswith('@SQ')]
    chrom_sizes = []
    chroms = []
    for i in sorted(chrom_lines, key=sort_chroms):
        for j in i.split('\t'):
            if j.startswith('LN'): chrom_sizes.append(int(j[3:]))
            if j.startswith('SN'): chroms.append(int(re.search('\d+', j[3:]).group()))
    return chroms, chrom_sizes

################
# variant call
################
def key_by_pos(x):
    gt, line = x
    arr = line.split('\t')
    chrm = int(re.search('\d+',arr[2]).group())
    pos = int(arr[3])
    return ((chrm, pos), (gt, line))

def no_unaln(x):
    gt, line = x
    if not line: return False
    arr = line.split('\t')
    return '*' not in arr[2]

def snp_call(aln_rdd, args, sam_files, ref, header, prefix):
    def get_bam(part_idx, records, samples, header):
        ret = []
        sam2fp={}
        for s in samples:
            sam_filename = "/mnt/%s.%d.sam" % (s, part_idx)
            ret.append(sam_filename)
            fw = open(sam_filename, 'w')
            sam2fp[s] = fw
            sam2fp[s].write(header)
        for rec in records:
            s, line = rec[1][1]
            sam2fp[s].write(line+'\n')
        for s in sam2fp: sam2fp[s].close()
        yield ','.join(ret)
    def call_snp(sam_files, ref):
        bam_files = []
        for sam in sam_files.split(','):
            bamfname = sam[:-4]+'.bam'
            bam_files.append(bamfname)
            cmd = "samtools view -bS %s | samtools sort > %s" % (sam, bamfname)
            call(cmd, shell=True)
            call("rm -f %s" % sam, shell=True)
        fid = (bam_files[0].split('.'))[-2]
        outfile = '/user/spark/warehouse/%s.%s.snp.txt' % (prefix, fid)
        snpcall_cmd = "/mnt/call_variants.py -b %s -f %s -u %d -c %d -p %f --hets_cov %d --hets_purity %f --indel_cov %d --indel_purity %f --indel_perc %f"\
        % (bam_files, ref, args['qual'], args['min_cov'], args['min_purity'], args['hets_cov'], args['hets_purity'], args['indel_cov'], args['indel_purity'], args['indel_perc'])
        if args['indel']: snpcall_cmd += " --indel"
        if args['hets']:  snpcall_cmd += " --hets"
        if args['vcf']:  snpcall_cmd += " --vcf"
        if args['mpileup_param']: snpcall_cmd += " -a '%s'" % args['mpileup_param']
        cmd = "%s | hadoop fs -put - %s" % (snpcall_cmd, outfile)
        call(cmd, shell=True)
        call("rm -f %s" % ' '.join(bam_files), shell=True)
        return outfile

    snp_files = aln_rdd.mapPartitionsWithIndex(lambda x, y: get_bam(x, y, sam_files, header))\
                       .map(lambda x: call_snp(x, ref))\
                       .collect()
    return snp_files

def snp_call_pipe(aln_rdd, args, sam_files, ref, header, prefix):
    def get_bam(part_idx, records, samples, header):
        ret = []
        bam2fp={}
        fw={}
        for s in samples:
            bam_filename = "/mnt/%s.%d.bam" % (s, part_idx)
            ret.append(bam_filename)
            fw[s] = open(bam_filename, 'w')
            cmd = "samtools view -bS - | samtools sort -"
            bam2fp[s] = Popen(cmd, stdin=PIPE, stdout=fw[s], bufsize=0, shell=True)
            #cmd = "samtools view -bS - | samtools sort - > %s" % bam_filename
            #bam2fp[s] = Popen(cmd, stdin=PIPE, bufsize=0, shell=True)
            bam2fp[s].stdin.write(bytes(header+'\n', 'UTF-8'))
        for rec in records:
            s, line = rec[1][1]
            bam2fp[s].stdin.write(bytes(line+'\n', 'UTF-8'))
        for s in bam2fp:
            bam2fp[s].communicate()
            bam2fp[s].stdin.close()
            fw[s].close()
        yield ','.join(ret)

    def call_snp(bam_files, ref):
        fid = ((bam_files.split(','))[0].split('.'))[-2]
        outfile = '/user/spark/warehouse/%s.%s.snp.txt' % (prefix, fid)
        snpcall_cmd = "/mnt/call_variants.py -b %s -f %s -u %d -c %d -p %f --hets_cov %d --hets_purity %f --indel_cov %d --indel_purity %f --indel_perc %f"\
        % (bam_files, ref, args['qual'], args['min_cov'], args['min_purity'], args['hets_cov'], args['hets_purity'], args['indel_cov'], args['indel_purity'], args['indel_perc'])
        if args['indel']: snpcall_cmd += " --indel"
        if args['hets']:  snpcall_cmd += " --hets"
        if args['vcf']:  snpcall_cmd += " --vcf"
        if args['mpileup_param']: snpcall_cmd += " -a '%s'" % args['mpileup_param']
        cmd = "%s | hadoop fs -put - %s" % (snpcall_cmd, outfile)
        call(cmd, shell=True)
        #call("rm -f %s" % ' '.join(bam_files), shell=True)
        return outfile
        #return snpcall_cmd

    snp_files = aln_rdd.mapPartitionsWithIndex(lambda x, y: get_bam(x, y, sam_files, header))\
                       .map(lambda x: call_snp(x, ref))\
                       .collect()
    return snp_files

def partition_rdd(sam_rdd, chroms, chrm_size, npartitions=400):
    def to_cumulative_chrom_size(x):
        ret = [0]
        cum=0
        for i in x:
            cum += i
            ret.append(cum)
        return ret

    def to_cumulative(chrom, pos):
        ind = chroms.index(chrom)
        return chrm_size_cumulative[ind]+pos
 
    def rangePartitioner(x):
        cum = chrm_size_cumulative
        chrm, pos = x
        if len(cum) == 1: cns_pos = int(pos)
        else:
            cns_pos = to_cumulative(chrm, pos)
        p = bisect.bisect_left(bounds, cns_pos)
        return p

    def get_bounds():
        nrecords_part = int(nrecords/npartitions)
        c = 0
        bounds_nonoverlap = []
        for chrm, pos in sam_rdd.sortByKey().map(lambda x: x[0]).collect():
            c+=1
            cns_pos = to_cumulative(chrm, pos)
            if c%nrecords_part == 0: bounds_nonoverlap.append(cns_pos)

        ret=[(1, bounds_nonoverlap[0])]
        for i in range(len(bounds_nonoverlap)-1):
            s = bounds_nonoverlap[i] - 300
            e = bounds_nonoverlap[i+1]
            ret.append((s, e))
        return ret

    def flatmap_bounds(bounds):
        if len(bounds) == 1: return [bounds[0][0], bounds[0][1]]
        coords=[]
        for s, e in bounds:
            coords.append(s)
            coords.append(e)
        coords.sort()
        return coords

    def toindexed(x, y):
        for d in y: yield (x, d)

    def myflatmap(x):
        ind, data = x
        if ind == 0: return [x]
        if ind%2 == 0:
            return [(ind-1, data), (ind+1, data)]
        else:
            return [x]

    def myPartitioner(x):
        return int(x)


    nrecords = sam_rdd.count()
    genome_size = sum(chrm_size)
    chrm_size_cumulative = to_cumulative_chrom_size(chrm_size)
    bounds = get_bounds()
    bounds = flatmap_bounds(bounds)
    return sam_rdd.partitionBy(npartitions*2, rangePartitioner)\
                  .mapPartitionsWithIndex(toindexed)\
                  .flatMap(myflatmap)\
                  .partitionBy(npartitions+1, myPartitioner)

def merge_SNP(args, snp_files, output_path, outfile, prefix):
    # Copy partial SNP files in HDFS to local folder
    call('hadoop fs -copyToLocal /user/spark/warehouse/* /mnt/', shell=True)
    call('hdfs dfs -rm /user/spark/warehouse/*', shell=True)

    # Create list of partial SNP files
    tmp_files = os.listdir('/mnt/')
    snps = ['/mnt/{}'.format(x) for x in tmp_files if x[-8:] == '.snp.txt']

    # output path
    outpath_local = '/mnt/'+outfile
    fw = open(outpath_local, 'w')

    if args['vcf']:
        header = check_output("grep '^#' %s" % ('/mnt/'+os.path.basename(snp_files[0])), shell=True).decode("utf-8")
        fw.write(header)
        number_header_lines = len(header.strip().split('\n'))
        merge_cmd = 'tail -n +%d -q %s | sort -k1V,1 -k2n,2 -k8V,8 | uniq' % (number_header_lines, ' '.join(snps))
        snp_line = Popen(merge_cmd, stdout=PIPE, shell=True)
        arr = snp_line.communicate()[0].decode("utf-8").strip().split('\n')
        nsnp = len(arr)
        for i in range(nsnp-1):
            this_arr = arr[i].split('\t')
            next_arr = arr[i+1].split('\t')
            this_pos = this_arr[1]
            next_pos = next_arr[1]
            this_ref = this_arr[3]
            next_ref = next_arr[3]
            this_alt = this_arr[4]
            next_alt = next_arr[4]
            if this_pos == next_pos and this_ref == next_ref and this_alt == next_alt: continue
            fw.write(arr[i]+'\n')
        if nsnp > 1: fw.write(arr[i+1]+'\n')
        fw.close()
    else:
        header = check_output('head -1 %s' % ('/mnt/'+os.path.basename(snp_files[0])), shell=True).decode("utf-8")
        header = header.replace('.0', '')
        fw.write(header)
        header = header.strip().split('\t')
        nsamples = int((len(header)-4)/3)
        sortStr = ''
        c=6
        for i in range(nsamples):
            sortStr = ' -k%sn,%s' % (str(c), str(c))
            c+=3
        merge_cmd = 'tail -n +2 -q %s | sort -k1V,1 -k2n,2%s | uniq' % (' '.join(snps), sortStr)
    
        snp_line = Popen(merge_cmd, stdout=PIPE, shell=True)
        arr = snp_line.communicate()[0].decode("utf-8").strip().split('\n')
        nsnp = len(arr)
        for i in range(nsnp-1):
            this_arr = arr[i].split('\t')
            next_arr = arr[i+1].split('\t')
            this_pos = this_arr[1]
            next_pos = next_arr[1]
            if this_pos == next_pos: continue
            fw.write(arr[i]+'\n')
        if nsnp > 1: fw.write(arr[i+1]+'\n')
        fw.close()

    # copy to s3 bucket
    s3_cp_cmd = "aws s3 cp %s %s" % (outpath_local, output_path+outfile)
    call(s3_cp_cmd, shell=True)

def save_BAMs(sam_rdd, header, output_path, gts, index):
    def by_sample(x, gt2index):
        return gt2index[x]

    def write_BAMs(recs, header, output_path):
        for gt, line in recs:
            first_line = line
            outfile = output_path+gt+'.bam'
            break
        cmd = "samtools view -bS - | samtools sort - | aws s3 cp - {}".format(outfile)
        bam = Popen(cmd, stdin=PIPE, bufsize=0, shell=True)
        bam.stdin.write(bytes(header+'\n', 'UTF-8'))
        bam.stdin.write(bytes(first_line+'\n', 'UTF-8'))
        for gt, line in recs:
            if not line: break
            bam.stdin.write(bytes(line+'\n', 'UTF-8'))
        bam.communicate()
        bam.stdin.close()
        yield outfile

    def write_indexed_BAMs(recs, header, output_path):
        for gt, line in recs:
            first_line = line
            outfile = output_path+gt+'.bam'
            outfile_local = '/mnt/'+gt+'.bam'
            break
        fw = open(outfile_local, 'w')
        cmd = "samtools view -bS - | samtools sort -"
        bam = Popen(cmd, stdin=PIPE, stdout=fw, bufsize=0, shell=True)
        bam.stdin.write(bytes(header+'\n', 'UTF-8'))
        bam.stdin.write(bytes(first_line+'\n', 'UTF-8'))
        for gt, line in recs:
            if not line: break
            bam.stdin.write(bytes(line+'\n', 'UTF-8'))
        bam.communicate()
        bam.stdin.close()
        fw.close()

        cmd = "samtools index {}; aws s3 cp {} {}; aws s3 cp {} {}".format(outfile_local, outfile_local, outfile, outfile_local+'.bai', outfile+'.bai')
        call(cmd, shell=True)
        yield outfile
   
 
    sample_size = len(gts)
    gt2index={}
    for i, v in enumerate(gts): gt2index[v] = i
    if index: 
        return sam_rdd.partitionBy(sample_size, lambda x: by_sample(x, gt2index))\
                      .mapPartitions(lambda x: write_indexed_BAMs(x, header, output_path))\
                      .collect()
    else:
        return sam_rdd.partitionBy(sample_size, lambda x: by_sample(x, gt2index))\
                      .mapPartitions(lambda x: write_BAMs(x, header, output_path))\
                      .collect()

def get_min_partitions(executor_memory, executor_cores):
    with open('/mnt/var/lib/info/job-flow.json', 'r') as f:
        j = json.load(f)
        instance_groups = j["instanceGroups"]
        for group in instance_groups:
            if group["instanceRole"] == "Core":
                nodes = group["requestedInstanceCount"]
                node_type = group["instanceType"]

    instance_types = {
        'm1.xlarge': {'vcpu': 4, 'memory': 12288},
        'cc2.8xlarge': {'vcpu': 32, 'memory': 56320},
        'c3.8xlarge': {'vcpu': 32, 'memory': 53248},
        'c4.8xlarge': {'vcpu': 36, 'memory': 53248},
        'm4.10xlarge': {'vcpu': 40, 'memory': 155648},
        'm4.16xlarge': {'vcpu': 64, 'memory': 253952},
    }
    instance_executors = int(math.floor((0.9 * instance_types[node_type]['memory']/1024) / executor_memory )) - 2
    min_partitions = instance_executors * nodes * executor_cores
    return min_partitions

# Main routine
def main(start_time):
    # Get command line options
    args = parse_args()
    avros = args['sample'].split(',')
    gts = [os.path.basename(fname)[:-5] for fname in avros]
    if args['vcf']:
        outfile = args['job_id']+'.vcf'
    else:
        outfile = args['job_id']+'.txt'
    output_path = args['output']
    prefix = args['job_id']

    # get header
    header=get_sam_header(args)

    # get chr sizes
    chroms, chrom_size = get_chrom_info(header)
    
    # Configure Spark cluster
    conf = SparkConf()\
           .setAppName("spark-aln")\
           .set("spark.executor.extraClassPath","/usr/lib/spark/examples/jars/spark-examples.jar")\
           .set("spark.serializer", "org.apache.spark.serializer.KryoSerializer")\
           .set("spark.akka.frameSize", 256)\
           .set("spark.python.worker.memory", "1g")
    sc = SparkContext(conf = conf)

    # Get job configuration
    executor_memory = int(re.findall('\d+', SparkConf().get("spark.executor.memory"))[0])
    executor_cores = int(SparkConf().get("spark.executor.cores"))

    # Make sure the AVRO file is partitioned into at least as many parts as the number of cluster executors
    min_partitions = get_min_partitions(executor_memory, executor_cores)

    # Get list of files to process
    input_samples = args['sample']
    samples = input_samples.split(',')

    tpartitions = 0
    sample2header={}
    sample_rdds = []
    for sample in samples:
        sam_rdd, nspartitions = reads2aln(sc, sample, args, min_partitions, executor_memory, executor_cores)
        sample_rdds.append(sam_rdd)
        tpartitions += nspartitions
    sam_rdd = sc.union(sample_rdds)
    del sample_rdds
   
    print('\n\n\tSamples are being processed as {:,} partitions\n\n'.format(tpartitions)) 
    # option to save intermediate BAM files
    if args['save_BAMs'] or args['aln_only']:
        BAM_files = save_BAMs(sam_rdd, header, output_path, gts, args['index'])
        if args['aln_only']: sys.exit()

    # keep only aligned records and use physical coordniates as key ((chrm, pos), (filename, record))
    sam_rdd = sam_rdd.filter(no_unaln)\
                     .map(key_by_pos)
    sam_rdd.cache()

    # range partition into chromosome segments
    sam_partition_rdd = partition_rdd(sam_rdd, chroms, chrom_size, npartitions=tpartitions)

    # genotype calling
    snp_files = snp_call_pipe(sam_partition_rdd, args, gts, '/mnt/'+args['ref'], header, prefix)

    print('\n\n\tMerging SNPs..\n\n')
    # merge SNPs and copy result to S3 bucket
    merge_SNP(args, snp_files, output_path, outfile, prefix)

 
if __name__ == "__main__":
    start_time = time.time()
    main(start_time)
    print("---total time {} seconds ---".format(time.time() - start_time))
