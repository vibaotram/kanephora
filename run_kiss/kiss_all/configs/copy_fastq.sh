fq=$(cat /shared/projects/elai_most/data/fastq/fastq_2019_clean/vn_accessions.txt)

for f in $fq
do
    f1=$(ls /shared/projects/elai_most/data/fastq/fastq_2019_clean/${f}/*1.fq.gz)
    ln -s ${f1} ${f}_R1.fastq.gz
    f2=$(ls /shared/projects/elai_most/data/fastq/fastq_2019_clean/${f}/*2.fq.gz)
    ln -s ${f2} ${f}_R2.fastq.gz
done