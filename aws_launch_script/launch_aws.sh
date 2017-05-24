#!/bin/bash

set -e

AWS_ACCESS_KEY=$1
AWS_SECRET_KEY=$2
CHROM=$3
PART=$4
KEYNAME=$5

usage()
{
    BASE=$(basename -- "$0")
    echo "Launch amazon instance to get diploid consensus sequence for each ssc bam parent
Usage:
    $BASE <aws_access_key> <aws_secret_key> <chrosome> <part> <keyname>
"
    exit 1
}
test -z ${AWS_ACCESS_KEY} && usage
test -z ${AWS_SECRET_KEY} && usage
test -z ${CHROM} && usage
test -z ${PART} && usage
test -z ${KEYNAME} && usage

# Instance details
SPOT_PRICE=0.05
INSTANCE_TYPE=c4.large
IMAGE_ID=ami-80861296

STARTUP_SCRIPT=$(cat /Users/shubhamsaini/Desktop/gatk-scripts/run_from_aws.sh | \
    sed "s/\$1/${AWS_ACCESS_KEY}/" | sed "s~\$2~${AWS_SECRET_KEY}~" | \
    sed "s~\$3~${CHROM}~" | \
    sed "s/\$4/${PART}/")
STARTUP_SCRIPT_ENCODE="$(echo "${STARTUP_SCRIPT}" | base64 -w 0)"

LAUNCH_SPEC="{\"ImageId\":\"${IMAGE_ID}\",\"Placement\":{\"AvailabilityZone\": \"us-east-1b\"},\"SecurityGroupIds\":[\"sg-5e914222\"], \"KeyName\":\"${KEYNAME}\",\"InstanceType\":\"${INSTANCE_TYPE}\", \"UserData\":\"${STARTUP_SCRIPT_ENCODE}\", \"BlockDeviceMappings\": [ {\"VirtualName\": \"string\",\"DeviceName\": \"/dev/sdf\",\"Ebs\": {\"VolumeSize\": 50,\"DeleteOnTermination\": true,\"VolumeType\": \"gp2\"}}]}"

aws ec2 request-spot-instances \
    --spot-price ${SPOT_PRICE} \
    --instance-count 1 \
    --type one-time \
    --launch-specification "${LAUNCH_SPEC}"