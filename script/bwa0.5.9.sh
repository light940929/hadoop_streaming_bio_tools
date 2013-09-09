#!/bin/bash
#coding :utf-8

Species ＝ Human_65.fa
getData = /home/hadoop/Data/getData
setData = /home/hadoop/Data/setData

$reference = http://192.168.1.198/hadoop/${Species}

wget -p $getData/${Species} $reference

$BAM_HOME/bwa index $getData/${Species}

$HADOOP_HOME/bin/hadoop fs -put ${Species}.amb  ${Species}.amb
$HADOOP_HOME/bin/hadoop fs -put ${Species}.ann  ${Species}.ann
$HADOOP_HOME/bin/hadoop fs -put ${Species}.bwt  ${Species}.bwt
$HADOOP_HOME/bin/hadoop fs -put ${Species}.pac  ${Species}.pac
$HADOOP_HOME/bin/hadoop fs -put ${Species}.rbwt ${Species}.rbwt
$HADOOP_HOME/bin/hadoop fs -put ${Species}.rpac ${Species}.rpac
$HADOOP_HOME/bin/hadoop fs -put ${Species}.rsa  ${Species}.rsa
$HADOOP_HOME/bin/hadoop fs -put ${Species}.sa   ${Species}.sa

rm -rf $getData/${Species}.*

$dna_pair_1 = testdata_1_1 
$dna_pair_2 = testdata_2_1

source = http://192.168.1.198/hadoop/${dna_pair}.fastq.gz

wget -p $setData/human_case   $source

for f in $( ls $setData/human_case); do
      eval gunzip "$setData/human_case/$f"
      $HADOOP_HOME/bin/hadoop fs -put $setData/human_case/$f $f
done

touch spcace.txt|$HADOOP_HOME/bin/hadoop fs -put spcace.txt spcace.txt|rm -rf space.txt


$HADOOP_HOME/bin/hadoop jar $HADOOP_JAR \
-files hdfs://192.168.1.198:54310/user/hadoop/$dna_pair_1.fq,hdfs://192.168.1.198:54310/user/hadoop/${Species},hdfs://192.168.1.198:54310/user/hadoop/${Species}.bwt,hdfs://192.168.1.198:54310/user/hadoop/${Species}.rbwt \
-input /user/hadoop/space.txt \
-output sai1.sai \
-mapper "$BAM_HOME/bwa aln dna_71_5000.fa testdata_1.fq"  \
-reducer NONE

$HADOOP_HOME/bin/hadoop jar $HADOOP_JAR \
-files hdfs://192.168.1.198:54310/user/hadoop/$dna_pair_2.fq,hdfs://192.168.1.198:54310/user/hadoop/${Species},hdfs://192.168.1.198:54310/user/hadoop/${Species}.bwt,hdfs://192.168.1.198:54310/user/hadoop/${Species}.rbwt \
-input /user/hadoop/space.txt \
-output sai2.sai \
-mapper "$BAM_HOME/bwa aln dna_71_5000.fa testdata_2.fq" \
-reducer NONE


$HADOOP_HOME/bin/hadoop jar $HADOOP_JAR \
-files hdfs://192.168.1.198:54310/user/hadoop/$dna_pair_1.fq,hdfs://192.168.1.198:54310/user/hadoop/$dna_pair_2.fq,hdfs://192.168.1.198:54310/user/hadoop/${Species},hdfs://192.168.1.198:54310/user/hadoop/sai1.sai,hdfs://192.168.1.198:54310/user/hadoop/sai2.sai,hdfs://192.168.1.198:54310/user/hadoop/${Species}.ann,hdfs://192.168.1.198:54310/user/hadoop/${Species}amb,hdfs://192.168.1.198:54310/user/hadoop/${Species}.pac,hdfs://192.168.1.198:54310/user/hadoop/${Species}.bwt,hdfs://192.168.1.198:54310/user/hadoop/${Species}.rbwt,hdfs://192.168.1.198:54310/user/hadoop/${Species}.sa,hdfs://192.168.1.198:54310/user/hadoop/${Species}.rsa \
-input /user/hadoop/space.txt \
-output aln.sam \
-mapper "$BAM_HOME/bwa sampe ${Species} sai1.sai sai2.sai  $dna_pair_1.fq $dna_pair_2.fq  " \
-reducer NONE