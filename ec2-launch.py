#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Launch AWS EC2 instances for converting fastq to avro.
"""

import sys, os
import boto3
import argparse
from subprocess import Popen, STDOUT

def parse_args():
    ap = argparse.ArgumentParser()
    ap.add_argument("sample_file", help="a file of list of samples (s3 object) to process (one sample for each line)")
    ap.add_argument("meta", help="a file of meta data")
    ap.add_argument("output", help="s3 object path for output")
    ap.add_argument("log", help="s3 object path for log")
    ap.add_argument("-n", metavar='INT', default=1, type=int, help="number of ec2 instances to create [1]")
    ap.add_argument("-b", metavar='STR', required=True, help="Bootstap script running on instance at launch [None]")
    args = ap.parse_args()
    return args

def create_ec2_instance(n, bootstrap):
    session = boto3.Session(profile_name='nonprod')
    ec2 = session.resource('ec2')
    #ec2 = session.client('ec2')

    if bootstrap == '':
        user_data = ''
    else:
        with open(bootstrap) as bs:
            user_data = bs.read()

    instance = ec2.create_instances(
        ImageId='ami-428aa838', # amazon linux 2
        #ImageId='ami-97785bed', # amazon base linux
        MinCount=n,
        MaxCount=n,
        InstanceType='t2.xlarge',
        KeyName='hui.guo',
        UserData=user_data,
        IamInstanceProfile={
            'Arn':  'arn:aws:iam::820225517789:instance-profile/DV-THES-Profile',
        },
        InstanceInitiatedShutdownBehavior='stop',
        NetworkInterfaces=[
            {   #'NetworkInterfaceId': 'eni-46b9cf9a',
                'SubnetId': 'subnet-f5d203af',
                'Groups': ['sg-16380368'],
                'AssociatePublicIpAddress': True,
                'DeleteOnTermination': True,
                'DeviceIndex': 0}
        ],
        BlockDeviceMappings=[
            {
                'DeviceName': '/dev/xvda',
                'Ebs': {
                    'DeleteOnTermination': True,
                    'SnapshotId': 'snap-0ef38d9386a70b23b', # amazon linux 2
                    #'SnapshotId': 'snap-0fae6f7252388fc12', # amazon base linux
                    'VolumeSize': 200,
                    'VolumeType': 'gp2'
                 }
            }
        ],
        TagSpecifications=[
            {
                'ResourceType': 'instance',
                'Tags': [
                    { 'Key': 'Name',    'Value': 'Fastq2avro' },
                ]
            }
        ]
    )
    return instance

def partition_samples(sample_file, n):
    with open(sample_file) as fp:
        samples = fp.read().splitlines()
    L = len(samples)
    assert 0 < n <= L
    s, r = divmod(L, n)
    t = s + 1
    return ([samples[p:p+t] for p in range(0, r*t, t)] +
            [samples[p:p+s] for p in range(r*t, L, s)])

def run_job(samples, instances, s3_path):
    session = boto3.Session(profile_name='nonprod')
    ec2 = session.client('ec2')
    ssm_client = session.client('ssm')
    # waiter
    waiter = ec2.get_waiter('instance_status_ok')
    # execute command
    for i, sample in zip(instances, samples):
        s = ','.join(sample)
        commands = ['fastq2avro_batch.py %s %s' % (s, s3_path)]
        instance_ids = [i.id]
        waiter.wait(InstanceIds=[i.id])
        resp = ssm_client.send_command(
            DocumentName="AWS-RunShellScript", 
            Parameters={'commands': commands},
            InstanceIds=instance_ids
        )
        print(resp)

def run_job2(samples, instances, s3_path, meta_file, log):
    session = boto3.Session(profile_name='nonprod')
    ec2 = session.client('ec2')
    waiter = ec2.get_waiter('instance_status_ok')
    FNULL = open(os.devnull, "w")
    for i, sample in zip(instances, samples):
        s = ','.join(sample)
        #commands = 'fastq2avro_batch.py %s %s' % (s, s3_path)
        commands = 'fastq2avro_batch_tar.py %s %s %s %s' % (s, meta_file, s3_path, log)
        cmd = "ssh -i /ngsprod/guoh1/huiguo_nonprod.pem -o \"StrictHostKeyChecking no\" ec2-user@%s \"%s\"" % (i.private_ip_address, commands) 
        waiter.wait(InstanceIds=[i.id])
        print ("Send command..")
        p = Popen(cmd, shell=True, stdout=FNULL, stderr=STDOUT)
        print ("Done.")

def main():
    args = parse_args()
    instances = create_ec2_instance(args.n, args.b)
    for i in instances:
        print (i.id)
    samples = partition_samples(args.sample_file, args.n)
    run_job2(samples, instances, args.output, args.meta, args.log)

if __name__ == "__main__":
    main()
