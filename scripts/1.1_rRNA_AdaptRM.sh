#!/bin/bash -l

echo "run script for rna-seq-analysis Cleaning and QC"

#./1.0_rRNA_AdaptRM.sh > 1.0_rRNA_AdaptRM.sh.log 2>&1
#1. Raw Data rRNA and AdaptRM Assessment

work_dir=/home/gala0002/proj/proj_dir/
mkdir -p ${work_dir}1.0_sort-trim/
out_dir1=${work_dir}1.0_sort-trim/

################################

#SID=$(echo $SID | sed 's=/[^/]*$==;s/\.$//')
#SID=$1
#mkdir ${out_dir1}${SID}

#----
SORTMERNADIR=/bioinfo/sortmerna-2.1b
TRIMMOMATIC=/bioinfo/Trimmomatic-0.39

scripts=/bioinfo/sortmerna-2.1b/scripts/

cd ${out_dir1}

#http://www.usadellab.org/cms/?page=trimmomatic
#for reference only (less sensitive for adapters)
#java -jar trimmomatic-0.35.jar PE -phred33 input_forward.fq.gz input_reverse.fq.gz output_forward_paired.fq.gz output_forward_unpaired.fq.gz output_reverse_paired.fq.gz output_reverse_unpaired.fq.gz ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

#P8403_303/02-FASTQ/180209_D00415_0305_ACC0NAANXX/
#for feven in 202 203 206

#do

#echo "SampleID:" ${feven}

for f in `ls /home/gala0002/proj/proj_dir/RNAseq_rawdata/NG-26198_F1*_1.fastq.gz`; do

DIR=$(dirname $f)
BASEFL=$(basename $f )
BASEFL_CP=$(basename $f | cut -d_ -f1,2,3,4,5,6 )
SID=$(basename $f | cut -d_ -f2,3,4,5,6 )

#$TRIMMOMATIC/adapters/
#NexteraPE-PE.fa, TruSeq3-PE-2.fa

echo "Dir_ID:" ${DIR}
echo "Sample_filepath:" ${BASEFL}
echo "Sample_ID:" ${BASEFL_CP}
echo "SID:" ${SID}

mkdir -p ${out_dir1}${SID}
temp_dir=${out_dir1}${SID}"/"


#3 Remove rRNA

gunzip -c ${DIR}"/"${BASEFL_CP}"_1.fastq.gz" > ${temp_dir}${SID}_1.fq
gunzip -c ${DIR}"/"${BASEFL_CP}"_2.fastq.gz" > ${temp_dir}${SID}_2.fq

${scripts}merge-paired-reads.sh ${temp_dir}${SID}_1.fq ${temp_dir}${SID}_2.fq ${temp_dir}${SID}-interleaved.fq

nice -n 5 ${SORTMERNADIR}/sortmerna --ref ${SORTMERNADIR}/rRNA_databases/silva-bac-16s-id90.fasta,${SORTMERNADIR}/index/silva-bac-16s-db:\
${SORTMERNADIR}/rRNA_databases/silva-bac-23s-id98.fasta,${SORTMERNADIR}/index/silva-bac-23s-db:\
${SORTMERNADIR}/rRNA_databases/silva-arc-16s-id95.fasta,${SORTMERNADIR}/index/silva-arc-16s-db:\
${SORTMERNADIR}/rRNA_databases/silva-arc-23s-id98.fasta,${SORTMERNADIR}/index/silva-arc-23s-db:\
${SORTMERNADIR}/rRNA_databases/silva-euk-18s-id95.fasta,${SORTMERNADIR}/index/silva-euk-18s-db:\
${SORTMERNADIR}/rRNA_databases/silva-euk-28s-id98.fasta,${SORTMERNADIR}/index/silva-euk-28s:\
${SORTMERNADIR}/rRNA_databases/rfam-5s-database-id98.fasta,${SORTMERNADIR}/index/rfam-5s-db:\
${SORTMERNADIR}/rRNA_databases/rfam-5.8s-database-id98.fasta,${SORTMERNADIR}/index/rfam-5.8s-db \
--reads ${temp_dir}${SID}-interleaved.fq --aligned ${temp_dir}${SID}-sort_aligned-rRNA \
--other ${temp_dir}${SID}-sort --log -v --paired_in -a 65 --fastx


# deinterleave a paired-end/matepair FASTQ file to the named pipes (IN BACKGROUND)
#bash /home/gala0002/proj/RNAseq-analysis-run/scripts_wheat_P15505_2020/deinterleave_fastq.sh P15505_102_sort_1.fq P15505_102_sort_2.fq < P15505_102_sort.fq
#${scripts}unmerge-paired-reads.sh ${temp_dir}${SID}-sort.fq ${temp_dir}${SID}-sort_1.fq ${temp_dir}${SID}-sort_2.fq

cat ${temp_dir}${SID}-sort.fq | grep '^@.* 1:' -A 3 --no-group-separator > ${temp_dir}${SID}-sort_1.fq
cat ${temp_dir}${SID}-sort.fq | grep '^@.* 2:' -A 3 --no-group-separator > ${temp_dir}${SID}-sort_2.fq

find . -name "${SID}-sort_[1,2].fq" | xargs -P 2 -I {} gzip {}


echo "Done zipping " ${temp_dir}${SID}

#4. Quality trimming and Adapter removal
#/usr/lib/jvm/java-8-openjdk-amd64/bin/java

nice -n 5 java -jar $TRIMMOMATIC/trimmomatic-0.39.jar PE -threads 65 -phred33 -trimlog ${temp_dir}trim_${SID}_logFile \
${temp_dir}${SID}-sort_1.fq.gz ${temp_dir}${SID}-sort_2.fq.gz \
${temp_dir}${SID}-sort-trim_1.fq.gz ${temp_dir}${SID}-sort-trim-unpaired_1.fq.gz \
${temp_dir}${SID}-sort-trim_2.fq.gz ${temp_dir}${SID}-sort-trim-unpaired_2.fq.gz \
ILLUMINACLIP:$TRIMMOMATIC/adapters/"TruSeq3-PE-2.fa":2:30:10 SLIDINGWINDOW:4:15 MINLEN:30


echo "Done sortmerna, trimommatic...."


rm ${temp_dir}${SID}_*.fq
rm ${temp_dir}${SID}-interleaved.fq
rm ${temp_dir}trim_${SID}_logFile
rm ${temp_dir}${SID}-sort.fq
rm ${temp_dir}${SID}-sort-trim-unpaired_*
rm ${temp_dir}${SID}-sort_*

echo "Done script sample > "${SID}

done

#done
########################################################################
echo "script all done...."

