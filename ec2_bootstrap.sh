#!/bin/bash

# make sure yum is up to date
sudo yum -y update

# install yum-utils, a collection of utilities and plugins that extend and supplement yum
sudo yum -y install yum-utils

# install development tools, which allow you to build and compile software from source code
sudo yum -y groupinstall development
sudo yum -y install openssl openssl-devel
sudo yum install -y gcc-c++ patch readline readline-devel zlib zlib-devel libyaml-devel libffi-devel openssl-devel make bzip2 autoconf automake libtool bison iconv-devel

# install python3 (amazon linux2)
sudo amazon-linux-extras install python3

# install boto3
sudo pip3.6 install boto3

# install java 8
sudo yum -y install java-1.8.0-openjdk.x86_64

# copy script to local
sudo aws s3 cp --recursive s3://dv-thes-ext/application/fastq2avro/ /usr/local/bin/
sudo chmod a+x /usr/local/bin/fastq2avro
sudo chmod a+x /usr/local/bin/fastq2avro_batch.py
sudo chmod a+x /usr/local/bin/fastq2avro_batch_tar.py

# copy meta data
aws s3 cp s3://dv-thes-ext/soy/soy_fq_csv.csv /home/ec2-user/

# install seqtk
git clone https://github.com/lh3/seqtk.git
cd seqtk; make
sudo cp seqtk /usr/local/bin/
cd ..
sudo rm -r seqtk

#### install SSM agent ####
# install go
#wget https://dl.google.com/go/go1.9.3.linux-amd64.tar.gz
#tar -zxf go1.9.3.linux-amd64.tar.gz
#sudo mv go/bin/* /usr/local/bin/
#rm -r *

# install rpm-build
#sudo yum install -y rpmdevtools rpm-build

# install ssm
#sudo yum install -y https://s3.amazonaws.com/ec2-downloads-windows/SSMAgent/latest/linux_amd64/amazon-ssm-agent.rpm
#sudo systemctl start amazon-ssm-agent
