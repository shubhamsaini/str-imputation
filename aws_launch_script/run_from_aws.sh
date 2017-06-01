#!/bin/bash

AWS_ACCESS_KEY=$1
AWS_SECRET_KEY=$2
CHROM=$3
PART_START=$4
PART_END=$5
BATCH_SIZE=$6


HOMEDIR=/root/
AWS_DIR=${HOMEDIR}/.aws
AWS_CONFIG_FILE=${AWS_DIR}/config
AWS_CRED_FILE=${AWS_DIR}/credentials


OUTBUCKET=s3://gatk-results

usage()
{
    BASE=$(basename -- "$0")
    echo "Run HipSTR
Usage:
    $BASE <aws access key> <aws secret key> <chrom> <part_start> <part_end> <batch_size>
       - aws access key and aws secret keys are for AWS configuration
       - chrom you are calling from
       - part_start is the first part to process
       - part_end is the last part to process
       - batch_size is the number of str to call in each run
Does the following:
1. Set up AWS configuration
2. Download necessary files
3. Create jobs
4. Run those jobs
5. Upload results to S3 bucket
6. Terminate
"
    terminate
    exit 1
}

terminate() {
    INSTANCE_ID=$(curl http://169.254.169.254/latest/meta-data/instance-id)
    # Get log
    aws s3 cp --output table /var/log/cloud-init-output.log ${OUTBUCKET}/log/${INSTANCE_ID}.log
    # Terminate instance
    echo "Terminating instance ${INSTANCE_ID}"
    aws ec2 terminate-instances --output table --instance-ids ${INSTANCE_ID}
    exit 1 # shouldn't happen
}

test -z ${AWS_ACCESS_KEY} && usage
test -z ${AWS_SECRET_KEY} && usage
test -z ${CHROM} && usage
test -z ${PART_START} && usage
test -z ${PART_END} && usage
test -z ${BATCH_SIZE} && usage


die()
{
    BASE=$(basename -- "$0")
    echo "$BASE: error: $@" >&2
    terminate
    exit 1
}

# Install things
sudo apt-get update || die "Could not update"
sudo apt-get -y install awscli || die "Could not install aws"
sudo apt-get -y install git || die "Could not install git"
sudo apt-get -y install make gcc libz-dev libncurses5-dev libbz2-dev liblzma-dev libcurl3-dev libssl-dev autoconf g++ python2.7 parallel || die "Could not install devtools"


cd ${HOMEDIR}
git clone https://github.com/samtools/htslib
cd htslib
git log --pretty=format:'%h' -n 1
autoheader
autoconf
./configure --enable-libcurl
make
sudo make install

cd ${HOMEDIR}
wget https://github.com/samtools/bcftools/releases/download/1.4.1/bcftools-1.4.1.tar.bz2
tar -xvjf bcftools-1.4.1.tar.bz2
cd bcftools-1.4.1
make
sudo make install

cd ${HOMEDIR}
git clone https://github.com/samtools/samtools
cd samtools
git log --pretty=format:'%h' -n 1
autoconf -Wno-syntax 
./configure
make
sudo make install

cd ${HOMEDIR}
git clone https://github.com/tfwillems/HipSTR.git
cd HipSTR
make

# Set up AWS credentials
echo "Setting up AWS credentials in ${AWS_DIR}"
mkdir -p ${AWS_DIR} || die "Could not create AWS dir"
echo "[default]" > ${AWS_CONFIG_FILE} || die "Could not write to ${AWS_CONFIG_FILE}"
echo "output = table" >> ${AWS_CONFIG_FILE} || die "Could not write to ${AWS_CONFIG_FILE}"
echo "region = us-east-1" >> ${AWS_CONFIG_FILE}  || die "Could not write to ${AWS_CONFIG_FILE}"
echo "aws_secret_access_key = ${AWS_SECRET_KEY}" >> ${AWS_CONFIG_FILE} || die "Could not write to ${AWS_CRED_FILE}"
echo "aws_access_key_id = ${AWS_ACCESS_KEY}" >> ${AWS_CONFIG_FILE} || die "Could not write to ${AWS_CRED_FILE}"
echo "[default]" > ${AWS_CRED_FILE} || die "Could not write to ${AWS_CRED_FILE}"
echo "aws_secret_access_key = ${AWS_SECRET_KEY}" >> ${AWS_CRED_FILE} || die "Could not write to ${AWS_CRED_FILE}"
echo "aws_access_key_id = ${AWS_ACCESS_KEY}" >> ${AWS_CRED_FILE} || die "Could not write to ${AWS_CRED_FILE}"

export AWS_ACCESS_KEY_ID=${AWS_ACCESS_KEY}
export AWS_SECRET_ACCESS_KEY=${AWS_SECRET_KEY}

# Set ulimit
echo "fs.file-max = 13107" | sudo tee -a /etc/sysctl.conf
sudo sysctl -p
echo "* soft     nproc          13107" | sudo tee -a /etc/security/limits.conf
echo "* hard     nproc          13107" | sudo tee -a /etc/security/limits.conf
echo "* soft     nofile         13107" | sudo tee -a /etc/security/limits.conf
echo "* hard     nofile         13107" | sudo tee -a /etc/security/limits.conf
echo "root soft     nproc          13107" | sudo tee -a /etc/security/limits.conf
echo "root hard     nproc          13107" | sudo tee -a /etc/security/limits.conf
echo "root soft     nofile         13107" | sudo tee -a /etc/security/limits.conf
echo "root hard     nofile         13107" | sudo tee -a /etc/security/limits.conf
echo "session required pam_limits.so" | sudo tee -a /etc/pam.d/common-session
ulimit -n 13107
ulimit -Hn 13107
echo "ulimit:"
ulimit -n


# setup ebs
sudo mkfs -t ext4 /dev/xvdf
sudo mkdir /storage
sudo mount /dev/xvdf /storage/
sudo chmod 777 /storage/

# Download files
cd /storage/ || die "Could not go to storage dir"

aws s3 cp ${OUTBUCKET}/phased/shapeit.chr$CHROM.vcf.gz . || die "Could not download phased SNP VCF"
tabix -p vcf shapeit.chr$CHROM.vcf.gz

mkdir indices
aws s3 cp s3://gatk-results/indices.tar.gz .
tar xvzf indices.tar.gz

mkdir fasta
cd fasta/
wget --retry-connrefused --waitretry=1 ftp://ftp.ncbi.nlm.nih.gov/1000genomes/ftp/technical/reference/human_g1k_v37.fasta.fai
wget --retry-connrefused --waitretry=1 ftp://ftp.ncbi.nlm.nih.gov/1000genomes/ftp/technical/reference/human_g1k_v37.fasta.gz
gunzip human_g1k_v37.fasta.gz
cd /storage/

git clone https://github.com/shubhamsaini/str-imputation.git || die "Could not clone github repo"

for CURR_PART in `seq ${PART_START} ${PART_END}`;
do
echo "Currently doing Part"
echo ${CURR_PART}

cd /storage/
cp -r str-imputation/hipstr_template/ hipstr_run_$CHROM\_${CURR_PART}
#cp indices/*.bai hipstr_run_$CHROM\_${CURR_PART}
cd hipstr_run_$CHROM\_${CURR_PART}

# Create jobs
let "nlines=${CURR_PART}*${BATCH_SIZE}"
head -n $nlines str_regions_bed/HipSTR.chr$CHROM.txt | tail -n ${BATCH_SIZE} > HipSTR_regions.txt
rstart=`awk '{print $2}' HipSTR_regions.txt | head -1`
rend=`awk '{print $3}' HipSTR_regions.txt | tail -1`
echo '#!/bin/sh' > samtoolsCommand.sh
awk -v part="${CURR_PART}" -v chrome="$CHROM" -v rstart="$rstart" -v rend="$rend" '{print "samtools view -b s3://sscwgs/"$1"/BAM/Sample_"$2"/analysis/"$2".final.bam "chrome":"rstart"-"rend" > /storage/hipstr_run_"chrome"\_"part"/bam/"$2".bam"}' famPhase1.txt >> samtoolsCommand.sh
chmod 777 samtoolsCommand.sh
mv samtoolsCommand.sh /storage/indices/
cd /storage/indices/
parallel < samtoolsCommand.sh
#rm *.final.bai

cd /storage/hipstr_run_$CHROM\_${CURR_PART}/bam/

# Find if any BAM files are corrupt, download if any
ls *.bam > files.list
cat files.list | xargs -I% -n1 samtools quickcheck -v % > redo.txt
echo "Redownloading following BAM files"
cat redo.txt
sed 's/\./\t/' redo.txt | awk '{print $1}' | xargs -I% -n1 grep % ../famPhase1.txt > ../famRedo.txt
cd ../
echo '#!/bin/sh' > samtoolsCommand.sh
awk -v part="${CURR_PART}" -v chrome="$CHROM" -v rstart="$rstart" -v rend="$rend" '{print "samtools view -b s3://sscwgs/"$1"/BAM/Sample_"$2"/analysis/"$2".final.bam "chrome":"rstart"-"rend" > /storage/hipstr_run_"chrome"\_"part"/bam/"$2".bam"}' famRedo.txt >> samtoolsCommand.sh
chmod 777 samtoolsCommand.sh
mv samtoolsCommand.sh /storage/indices/
cd /storage/indices/
parallel < samtoolsCommand.sh
#rm *.final.bai

cd /storage/hipstr_run_$CHROM\_${CURR_PART}/bam/
./sort.index.sh
ls *.bam > files.list


# Run jobs
~/HipSTR/HipSTR --bam-files files.list --fasta /storage/fasta/human_g1k_v37.fasta --regions ../HipSTR_regions.txt --str-vcf ../hipstr_calls_$CHROM\_${CURR_PART}.vcf.gz --snp-vcf /storage/shapeit.chr$CHROM.vcf.gz --log ../hipstr_calls_$CHROM\_${CURR_PART}.log --stutter-out ../hipstr_calls_$CHROM\_${CURR_PART}.stutter
cd ../
aws s3 cp hipstr_calls_$CHROM\_${CURR_PART}.vcf.gz ${OUTBUCKET}/hipstr/
aws s3 cp hipstr_calls_$CHROM\_${CURR_PART}.log ${OUTBUCKET}/hipstr/
aws s3 cp hipstr_calls_$CHROM\_${CURR_PART}.stutter ${OUTBUCKET}/hipstr/

cd /storage/
rm -rf hipstr_run_$CHROM\_${CURR_PART}
done

terminate