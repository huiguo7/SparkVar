# SparkVar
SparkVar is designed for swiftly processing large amounts of DNA sequencing data for read alignment and variant calling. SparkVar creates an EMR cluster in AWS to process sequencing data stored in S3.

### Dependencies:
1. An AWS account (https://aws.amazon.com)
2. boto3 installed on local machine (https://boto3.readthedocs.io/en/latest)
3. python3 (https://www.python.org/downloads)

### Convert fastq to avro
$./fastq2avro --help  
Usage: fastq2avro [options]  
  
  -1, --read1 <read1-file>    A FASTQ file containing read1 either gzipped or uncompressed  
  -2, --read2 <read2-file>    A FASTQ file containing read2 either gzipped or uncompressed  
  --help                      prints usage  
  
Examples...  
$fastq2avro --read1 r1.fastq.gz > out.avro  
$fastq2avro --read1 r1.fastq.gz --read2 r2.fastq.gz > out.avro  


### How to run SparkVar:
1. Copy the "application" folder under a S3 object path, specified with "--pipeline" option.
   SparkVar and its dependencies will be installed automatically during bootstrap stage when launching EMR.    

2. Option to convert reads in fastq format to avro.
   SparkVar takes avro files as input (see avro schema). If your reads data are stored in fastq format, use "fastq2avro" to convert them to avro. You may use ec2-launch.py to create a number of ec2 instances for this step. 
   
3. Launch emr cluster (emr-launch.py) to run SparkVar jobs.
   Run "emr-launch.py -h" for a number of input options. Make sure to include "/" at the end of each specified s3 object path. An example comamnd is included in the "example" folder.


### Avro schema:  
{  
    “name”: “FASTQ”,  
    “type”: “record”,  
    “fields”: [  
         {“name”: “header”, “type”: “string”, “default”: “”},  
         {“name”: “seq1”,   “type”: “string”, “default”: “”},  
         {“name”: “qual1”,  “type”: “string”, “default”: “”},  
         {“name”: “seq2”,   “type”: “string”, “default”: “”},  
         {“name”: “qual2”,  “type”: “string”, “default”: “”},  
    ]  
}  
