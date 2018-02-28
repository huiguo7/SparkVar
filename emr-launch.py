#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Launch AWS EMR cluster for running SparkVar jobs.
"""

import sys
import os
import boto3
import botocore
import argparse
import json
import datetime
import math

# Parse arguments
def parse_args():
    Defaults = {'MinMapQual':30, 'MinCov':3, 'MinPurity':0.98, 'HetsCov':10, 'HetsPurity':0.7, \
                'IndelCov':1, 'IndelPurity':0.5, 'IndelPerc':0.1 }

    core_node_list = ['cc2.8xlarge', 'c3.4xlarge', 'c3.8xlarge', 'c4.4xlarge',
                      'c4.8xlarge', 'm4.10xlarge', 'm4.16xlarge', 'm1.xlarge']

    ap = argparse.ArgumentParser()

    # Common options
    group0 = ap.add_argument_group('Required')
    group0.add_argument("--job_id", metavar='STR', default='SparkVar_job',
        help="optional job id")
    group0.add_argument("--samples", metavar='STR', required=True,
        help="a file of list of samples (s3 object) to process (one sample for each line)")
    group0.add_argument("--reference", metavar='STR', required=True,
        help="name of reference sequence file (fasta)")
    group0.add_argument("--ref_s3_path", metavar='STR', required=True,
        help="s3 object path of reference sequence")
    group0.add_argument("--bwt_prefix", metavar='STR', required=True,
        help="Prefix of bowtie2 reference index")
    group0.add_argument("--bwt_s3_path", metavar='STR', required=True,
        help="s3 object path for bowtie2 reference index")
    group0.add_argument("--output", metavar='STR', required=True,
        help="s3 object path for output")
    group0.add_argument("--pipeline", metavar='STR', required=True,
        help="s3 object path for data processing pipeline")
    group0.add_argument("--ec2_subnet_id", metavar='STR', required=True,
        help="EC2 subnet ID for EMR cluster")
    group0.add_argument("--emr_managed_master_security_group", metavar='STR', required=True,
        help="The identifier of the Amazon EC2 security group for the master node")
    group0.add_argument("--emr_managed_slave_security_group", metavar='STR', required=True,
        help="The identifier of the Amazon EC2 security group for the slave nodes")
    group0.add_argument("--service_access_security_group", metavar='STR', required=True,
        help="The identifier of the Amazon EC2 security group for the Amazon EMR service to access clusters in VPC private subnets")

    # alignment options
    group1 = ap.add_argument_group('Alignment')
    #group1.add_argument("--lib_type", metavar='STR', default='paired_end', choices=['single_end', 'paired_end'],
    #    help="Type of DNA library, choices are [single_end, paired_end]. [paired_end]")
    group1.add_argument("--aln_only", default=False, action="store_true",
        help="Option to only perform alignment (force '--save_BAMs'). [False]")
    group1.add_argument("--bt2_param", metavar='\'STR\'',
        help="option to bowtie2 parameters.")

    # SNP-calling options
    group2 = ap.add_argument_group('Variant calling')
    group2.add_argument("--qual", metavar='INT', type=int, default=Defaults['MinMapQual'],
        help="Option to set the minimum mapping quality for SNP call. [{}]".format(Defaults['MinMapQual']))
    group2.add_argument("--indel", default=False, action='store_true',
        help="Flag to call INDELs too, defaults to not call indels. [False]")
    group2.add_argument("--hets", default=False, action='store_true',
        help="Flag to call heterozygous SNPs. [False]")
    group2.add_argument("--vcf", default=False, action='store_true',
        help="Option to call variants using samtools/bcftools. [tabular]")
    group2.add_argument("--mpileup_param", metavar='\'STR\'',
        help="Option to set samtools mpileup parameters.")
    group2.add_argument("--min_cov", metavar='INT', default=Defaults['MinCov'], type=int,
        help="Option to set the minimum coverage for snp call. [{}]".format(Defaults['MinCov']))
    group2.add_argument("--min_purity", metavar='FLOAT', default=Defaults['MinPurity'], type=float,
        help="Option to set the minimum purity for snp call. [{}]".format(Defaults['MinPurity']))
    group2.add_argument("--hets_cov", metavar='INT', default=Defaults['HetsCov'], type=int,
        help="Option to set the minimum coverage for hets call. [{}]".format(Defaults['HetsCov']))
    group2.add_argument("--hets_purity", metavar='FLOAT', default=Defaults['HetsPurity'], type=float,
        help="Option to set the maximum purity for hets call. [{}]".format(Defaults['HetsPurity']))
    group2.add_argument("--indel_cov", metavar='INT', default=Defaults['IndelCov'], type=int,
        help="Option to set the minimum coverage for indel call. [{}]".format(Defaults['IndelCov']))
    group2.add_argument("--indel_purity", metavar='FLOAT', default=Defaults['IndelPurity'], type=float,
        help="Option to set the minimum purity for indel call. [{}]".format(Defaults['IndelPurity']))
    group2.add_argument("--indel_perc", metavar='FLOAT', default=Defaults['IndelPerc'], type=float,
        help="Option to set the minimum percent for displaying an indel allele. [{}]".format(Defaults['IndelPerc']))

    # EMR Cluster configuration options
    group3 = ap.add_argument_group('EMR configuration:')
    group3.add_argument("--EMR_release", metavar='STR', default="emr-5.5.0",
        help="The release label for the Amazon EMR release. [emr-5.5.0]")
    group3.add_argument("--master_node_type", metavar='STR', default="m1.xlarge",
        help="EC2 instance type of the master node. [m1.xlarge]")
    group3.add_argument("--master_node_number", metavar='INT', default=1,
        help="Number of master nodes to assign to the EMR cluster. [1]")
    group3.add_argument("--master_node_ebs", metavar='INT', type=int,
        help="EBS volume attched to master nodes (Gb). EBS volume of 500 Gb is attached for EBS-only instance types [None]")
    group3.add_argument("--additional_master_security_groups", metavar='STR', default="",
        help="Additional Amazon EC2 security group IDs for the master node (delimited by ',' if multiple IDs were given). [None]")
    group3.add_argument("--core_node_type", metavar='STR', default="m4.10xlarge",
        help="EC2 instance type of the core nodes. Options are "+', '.join(core_node_list)+'. [m4.10xlarge]')
    group3.add_argument("--core_node_number", metavar='INT', default=8,
        help="Number of core nodes to assign to the EMR cluster. [8]")
    group3.add_argument("--core_node_ebs", metavar='INT', type=int,
        help="EBS volume attched to core nodes (Gb). EBS volume of 500 Gb is attached for EBS-only instance types [None]")
    group3.add_argument("--additional_slave_security_groups", metavar='STR', default="",
        help="Additional Amazon EC2 security group IDs for the slave nodes (delimited by ',' if multiple IDs were given). [None]")
    group3.add_argument("--ssh_key", metavar='STR', default=None,
        help="SSH key to use with master node. [None]")
    group3.add_argument("--emr_cluster_name", metavar='STR', default="My Cluster",
        help="Name used to identify the EMR cluster. [My Cluster]")
    group3.add_argument("--job_flow_role", metavar='STR', default="EMR_EC2_DefaultRole",
        help=" An IAM role for an EMR cluster. [EMR_EC2_DefaultRole]")
    group3.add_argument("--tag", metavar='STR', default="",
        help="Name used to tag EC2 instances assocaited with EMR cluster. [None]")
    group3.add_argument("--auto_terminate", default=False, action="store_true",
        help="Option to terminate the EMR cluster upon completion of all steps. [False]")
    group3.add_argument("--deploy_mode", metavar='STR', default='cluster', choices=['client', 'cluster'],
        help="Where to run the Spark application from. Options are [client, cluster]. [cluster]")
    group3.add_argument("--log", metavar='STR', default=None,
        help="s3 object path for log infomation. [None]")

    # Spark options
    group4 = ap.add_argument_group('Spark:')
    group4.add_argument("--executor_memory", metavar='STR',
        help="Amount of memory per executor (e.g. 5G). [inferred from type of instances]")
    group4.add_argument("--executor_cores", metavar='INT', type=int,
        help="Number of cores per executor. [inferred from type of instances]")
    group4.add_argument("--number_partitions", metavar='INT', type=int,
        help="Number of data partitions. [inferred from type of instances]")

    # File Storage options
    group5 = ap.add_argument_group('File storage:')
    group5.add_argument("--save_BAMs", default=False, action="store_true",
        help="Flag to save BAM files after processing is complete. [False]")
    group5.add_argument("--index", default=False, action="store_true",
        help="Option to index the final bam file. [False]")

    args = vars(ap.parse_args())

    if args['core_node_type'].lower() not in core_node_list:
        print ("Here are the choices for 'core_node_type':")
        for i in core_node_list: print (i)
        sys.exit(1)

    return args

# Function to build list of EMR applications to install
def build_applications():
    return [
        {
            'Name': 'Spark',
        },
        {
            'Name': 'Ganglia',
        },
    ]

# Function to build Bootstrap Actions
def build_bootstrap(args):
    # Build the Bootstrap JSON
    bootstrap_action = [
        {
            'Name': "Install Tools",
            'ScriptBootstrapAction': {
                'Path': "{}application/SparkVar/install_tools_bootstrap_param.sh".format(args['pipeline']),
                'Args': [args["pipeline"], args["ref_s3_path"], args['bwt_s3_path']]
            }
        }
    ]

    return bootstrap_action

# Function to build the configuration JSON for the EMR cluster
def build_configuration():
    return [
        {
            "Classification":"spark-env",
            "Properties":{},
            "Configurations":[
                {
                    "Classification":"export",
                    "Properties":{
                        "PYSPARK_PYTHON":"/usr/bin/python3",
                        "PYTHONHASHSEED":"0"
                    },
                    "Configurations":[]
                }
            ]
        },
        {
            "Classification":"spark",
            "Properties":{
                "maximizeResourceAllocation":"true"
            },
            "Configurations":[]
        },
        {
            "Classification":"hadoop-env",
            "Properties":{},
            "Configurations":[
                {
                    "Classification":"export",
                    "Properties":{
                        "PYTHONHASHSEED":"0"
                    },
                    "Configurations":[]
                }
            ]
        },
        {
            "Classification":"spark-defaults",
            "Properties":{
                "spark.yarn.appMasterEnv.PYTHONHASHSEED":"0"
            },
            "Configurations":[]
        }
    ]

# Function to build EMR EC2 Instance-Groups
def instance_groups(args):
    master = {
        'Name': 'Master - 1',
        'Market': 'ON_DEMAND',
        'InstanceRole': 'MASTER',
        'InstanceType': args["master_node_type"],
        'InstanceCount': int(float(args["master_node_number"])),
    }
    if args['master_node_ebs']:
        master['EbsConfiguration'] = ebs_config(args['master_node_ebs'])
    elif args['master_node_type'] in ['c4.4xlarge', 'c4.8xlarge', 'm4.10xlarge', 'm4.16xlarge']:
        master['EbsConfiguration'] = ebs_config(500)

    core = {
        'Name': 'Core - 2',
        'Market': 'ON_DEMAND',
        'InstanceRole': 'CORE',
        'InstanceType': args["core_node_type"],
        'InstanceCount': int(float(args["core_node_number"])),
    }
    if args['core_node_ebs']:
        core['EbsConfiguration'] = ebs_config(args['core_node_ebs'])
    elif args['core_node_type'] in ['c4.4xlarge', 'c4.8xlarge', 'm4.10xlarge', 'm4.16xlarge']:
        core['EbsConfiguration'] = ebs_config(500)

    return [master, core]

def ebs_config(volume_size):
    return  {
                'EbsBlockDeviceConfigs': [
                    {
                        'VolumeSpecification': {
                            'VolumeType': 'gp2',
                            #'Iops': 10000,
                            'SizeInGB': volume_size
                        },
                         'VolumesPerInstance': 1
                    },
                 ],
                 'EbsOptimized': False
            }
    
# Function to build Tags
def build_tags(args):
    if args["tag"] is None:
        return None
    else:
        return [
            {
                'Key': 'Name',
                'Value': args["tag"]
            }
        ]

# Function to set the Spark-Submit parameters based on the EC2 instance type and reference assembly
def get_spark_parameters(args):
    # Get the total size of BT2 files associated with the reference assembly
    session = boto3.Session(profile_name='nonprod')
    client = session.client('s3')
    
    r = client.list_objects_v2(Bucket=(args['bwt_s3_path'].split('/'))[2],
                               Prefix='/'.join((args['bwt_s3_path'].split('/'))[3:])
                               )

    # Sum up the sizes of the individual bt2 files
    size = 0
    for item in r['Contents']:
        if item['Key'].endswith('bt2'):
            size += item['Size']

    # Convert to GB
    bt2_size = round(float(size) / (1024**3),2)

    # List of memory availabe to EMR EC2 instances within a Spark application
    # (http://docs.aws.amazon.com/emr/latest/ReleaseGuide/emr-hadoop-task-config.html)
    instance_types = {
        'm1.xlarge': {'vcpu': 4, 'memory': 12288},
        'cc2.8xlarge': {'vcpu': 32, 'memory': 56320},
        'c3.8xlarge': {'vcpu': 32, 'memory': 53248},
        'c4.8xlarge': {'vcpu': 36, 'memory': 53248},
        'm4.10xlarge': {'vcpu': 40, 'memory': 155648},
        'm4.16xlarge': {'vcpu': 64, 'memory': 253952},
    }

    # Get memory and vCPU count for the instance type
    instance = args['core_node_type']
    available_memory = (float(instance_types[instance]['memory']) / 1024) / 1.1
    vcpus = instance_types[instance]['vcpu']
    blocking_memory = int(math.ceil(0.75 * instance_types[instance]['memory'] / 1024))

    # Find the number of executors to minimize the number of wasted vCPUs
    config = {'wasted_cpus': 1000, 'vcpus_per_executor': 0, 'executors': 0}
    for executor_memory in range(int(math.ceil(bt2_size)), int(math.floor(available_memory))):
        # Calculate the number of executors
        executors = math.floor(available_memory / (executor_memory))

        # Calculate how many cores are available to each executor
        cores = math.floor(vcpus / executors)

        # Calculate how many Bowtie2 instances can run on the executor
        bowties = math.floor(executor_memory / bt2_size)
        if bowties > 5:
            bowties = 5

        # Calculate how many Bowtie2 instances can run on the instance
        total_bowties = bowties * executors

        # Calculate the number of threads for each Bowtie instance
        threads = math.floor(cores / bowties)

        # Get the number of wasted vCPUs
        waste = vcpus - (threads * bowties * executors)

        # Update config if this configuration minimizes wasted vCPU and maximizes cores/executor
        if (waste < config['wasted_cpus']) and (cores > config['vcpus_per_executor']):
            excess_memory = available_memory - (executors * executor_memory)

            additional_per_executor = math.floor(excess_memory / executors)
            #excess_memory -= additional_per_executor

            config = {
                'wasted_cpus': waste,
                'vcpus_per_executor': cores,
                'executors': executors,
                'memory_per_executor': executor_memory + additional_per_executor,
                'bowties_per_executor': bowties,
                'vcpus_per_bowtie': threads,
                'bowtie_instances': total_bowties,
                'unused_memory': round(excess_memory,1)
            }

    # Return the selected --executor-memory and --executor-cores settings, and the memory setting required
    # for one executor to fit on the cluster
    return int(config['memory_per_executor']), int(config['bowties_per_executor']), blocking_memory

# Function to build Steps
def build_steps(args):
    executor_memory, executor_cores, blocking_memory = get_spark_parameters(args)

    if args['deploy_mode'] == 'cluster':
        SparkVar_path = '{}application/SparkVar/SparkVar.py'.format(args['pipeline'])
    else:
        SparkVar_path = '/mnt/SparkVar.py'

    # Get the list of samples
    with open(args['samples']) as fp:
        sample_list = fp.read().splitlines()
    samples = ','.join(sample_list)

    SparkVar_cmd = ["spark-submit",
                    "--deploy-mode", "{}".format(args['deploy_mode']),
                    "--master","yarn",
                    "--executor-memory","{}G".format(executor_memory),
                    "--executor-cores","{}".format(executor_cores),
                    "--jars","/usr/lib/spark/examples/jars/spark-examples.jar",
                    "{}".format(SparkVar_path),
                    "--job_id","{}".format(args['job_id']),
                    "--sample","{}".format(samples),
                    "--ref","{}".format(args["reference"]),
                    "--bt2","{}".format(args["bwt_prefix"]),
                    "--output","{}".format(args["output"]),
                    #"--lib_type","{}".format(args["lib_type"]),
                    "--qual","{}".format(args["qual"]),
                    "--min_cov","{}".format(args["min_cov"]),
                    "--min_purity","{}".format(args["min_purity"]),
                    "--hets_cov","{}".format(args["hets_cov"]),
                    "--hets_purity","{}".format(args["hets_purity"]),
                    "--indel_cov","{}".format(args["indel_cov"]),
                    "--indel_purity","{}".format(args["indel_purity"]),
                    "--indel_perc","{}".format(args["indel_perc"]),
                   ]
    if args['aln_only']: SparkVar_cmd += ["--aln_only"]
    if args['indel']: SparkVar_cmd += ["--indel"]
    if args['hets']: SparkVar_cmd += ["--hets"]
    if args['vcf']: SparkVar_cmd += ["--vcf"]
    if args['save_BAMs']: SparkVar_cmd += ["--save_BAMs"]
    if args['index']: SparkVar_cmd += ["--index"]
    if args['number_partitions']: SparkVar_cmd += ["--number_partitions", args['number_partitions']]
    if args['bt2_param']: SparkVar_cmd += ["--bt2_param", "'"+args['bt2_param']+"'"]
    if args['mpileup_param']: SparkVar_cmd += ["--mpileup_param", "'"+args['mpileup_param']+"'"]
    if args['executor_memory']: SparkVar_cmd[6:7] = [args['executor_memory']]
    if args['executor_cores']: SparkVar_cmd[8:9] = [str(args['executor_cores'])]
    
    steps = [
        {
            "Name": "SparkVar",
            "ActionOnFailure": "CONTINUE",
            'HadoopJarStep': {
                "Jar":"command-runner.jar",
                "Args": SparkVar_cmd
            }
        }
    ] 
    
    return steps

def main():
    # Get command line arguments
    args = parse_args()
    print (' '.join(build_steps(args)[0]['HadoopJarStep']['Args']))
        
    # Build EMR cluster and initiate steps
    session = boto3.Session(profile_name='nonprod')
    client = session.client('emr')

    kwargs = {'Applications' : build_applications(),
              'BootstrapActions': build_bootstrap(args),
              'Configurations' : build_configuration(),
              'Instances' : {
                               'InstanceGroups': instance_groups(args),
                               'Ec2KeyName': args["ssh_key"],
                               'KeepJobFlowAliveWhenNoSteps': (True if args["auto_terminate"]==False else False),
                               'TerminationProtected': False,
                               'Ec2SubnetId': args['ec2_subnet_id'],
                               'EmrManagedMasterSecurityGroup': args['emr_managed_master_security_group'],
                               'EmrManagedSlaveSecurityGroup': args['emr_managed_slave_security_group'],
                               'ServiceAccessSecurityGroup': args['service_access_security_group'],
                               'AdditionalMasterSecurityGroups': args['additional_master_security_groups'].split(','),
                               'AdditionalSlaveSecurityGroups': args['additional_slave_security_groups'].split(','),
                             },
              'JobFlowRole' : args['job_flow_role'],
              'Name' : args['emr_cluster_name'],
              'ReleaseLabel' : args['EMR_release'],
              'ScaleDownBehavior' : 'TERMINATE_AT_INSTANCE_HOUR',
              'ServiceRole' : args['job_flow_role'],
              'Steps' : build_steps(args),
              'Tags' : build_tags(args),
              'VisibleToAllUsers' : True
             }
    if args['log']: kwargs['LogUri'] = args['log']

    response = client.run_job_flow(**kwargs)

    if 'JobFlowId' in response:
        cluster_id = response['JobFlowId']

        # Print status message for user
        msg = "A cluster with {} x {} master node and {} x {} core nodes has been started. ".format(args['master_node_number'], args['master_node_type'],
                                                                                                   args['core_node_number'], args['core_node_type'])
        msg += "Your cluster id is {}.".format(cluster_id)

        print(msg)
    else:
        print("Your cluster failed to start. Sorry!")

# Main executable
if __name__ == "__main__":
    main()
