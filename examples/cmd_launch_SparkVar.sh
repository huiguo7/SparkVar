#!/bin/bash

python3 ../emr-launch.py --job_id SparkVar_test \
                 --samples input_samples.txt \
                 --reference Zea_mays.AGPv4.fa \
                 --ref_s3_path s3://dv-thes-ext/maize_ref/ \
                 --bwt_prefix Zea_mays.AGPv4.fa \
                 --bwt_s3_path s3://dv-thes-ext/maize_ref/ \
                 --output s3://dv-thes-ext/test_out/ \
                 --pipeline s3://dv-thes-ext/ \
                 --emr_cluster_name dv-thes-ext \
                 --ec2_subnet_id *** \
                 --emr_managed_master_security_group *** \
                 --emr_managed_slave_security_group *** \
                 --service_access_security_group *** \
                 --additional_master_security_groups *** \
                 --additional_slave_security_groups *** \
                 --job_flow_role *** \
                 --ssh_key hui.guo \
                 --core_node_number 5 \
                 --log s3://dv-thes-ext/log/ \
                 --save_BAMs \
                 --index

