# !/bin/bash

# Change directories
cd /mnt/

# make sure yum is up to date
sudo yum -y update

# install yum-utils, a collection of utilities and plugins that extend and supplement yum
sudo yum -y install yum-utils

# install development tools, which allow you to build and compile software from source code
sudo yum -y groupinstall development
sudo yum -y install openssl openssl-devel
sudo yum install -y gcc-c++ patch readline readline-devel zlib zlib-devel libyaml-devel libffi-devel openssl-devel make bzip2 autoconf automake libtool bison iconv-devel bzip2-devel-1.0.6-13.amzn2.x86_64 xz-devel-5.2.2-1.amzn2.x86_64

# install python3 (amazon linux2)
sudo amazon-linux-extras install python3

# install bowtie2
wget https://iweb.dl.sourceforge.net/project/bowtie-bio/bowtie2/2.3.4/bowtie2-2.3.4-linux-x86_64.zip
unzip bowtie2-2.3.4-linux-x86_64.zip
sudo mv bowtie2-2.3.4-linux-x86_64/bowtie2* /usr/local/bin/
sudo rm -r bowtie2-2.3.4-linux-x86_64*

# install htslib
git clone https://github.com/samtools/htslib.git
cd htslib
autoheader
autoconf
./configure
make
sudo make install
cd ..
sudo rm -r htslib

# add lib path
sudo sh -c "echo '/usr/local/lib' >> /etc/ld.so.conf"

# install samtools
git clone https://github.com/samtools/samtools.git
cd samtools
autoheader
autoconf -Wno-syntax
./configure
make
sudo make install
cd ..
sudo rm -r samtools

# install bcftools
wget https://github.com/samtools/bcftools/releases/download/1.6/bcftools-1.6.tar.bz2
tar jxf bcftools-1.6.tar.bz2
cd bcftools-1.6
./configure
make
sudo make install
cd ..
sudo rm -r bcftools-1.6*

# Get s3 object path of application and reference
app_path=$1
reference=$2
bwt=$3

# copy script to local
sudo aws s3 cp ${app_path}application/SparkVar/ /usr/local/bin/ --recursive --exclude "*.sh"
sudo chmod a+x /usr/local/bin/call_variants.py /usr/local/bin/fasta_extract
sudo chmod 777 /mnt 

# Download FASTA file
aws s3 cp $reference /mnt/ --recursive

# Download bowtie2 index files
aws s3 cp $bwt /mnt/ --recursive
