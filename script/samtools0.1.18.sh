#!/bin/bash
#coding :utf-8

Species ＝ Human_65.fa
getData = /home/hadoop/Data/getData

$dna_pair = testdata

$outputsam = ${dna_pair}.sam

$outputmergesam = ${dna_pair}_m.sam

$outputbam = ${dna_pair}.bam

$outputmergebam = ${dna_pair}_m.bam

$outputrmbam = ${dna_pair}_rmdup.bam

$HADOOP_HOME/bin/hadoop fs -getmerge $outputsam $outputmergesam && $HADOOP_HOME/bin/hadoop fs -put $outputmergesam $outputmergesam | rm -rf $outputmergesam

$SAMTOOLS_HOME/samtools faidx  ${Species}
$HADOOP_HOME/bin/hadoop fs -put ${Species}.fai ${Species}.fai

rm -rf $getData/${Species}.fai

$HADOOP_HOME/bin/hadoop jar $HADOOP_JAR \
-files hdfs://192.168.1.198:54310/user/hadoop/${Species},hdfs://192.168.1.198:54310/user/hadoop/${Species}.fai,hdfs://192.168.1.198:54310/user/hadoop/$outputmergesam \
-input /user/hadoop/space.txt \
-output /user/hadoop/$outputbam \
-mapper "$SAMTOOLS_HOME/samtools view -bt ${Species}.fai $outputmergesam" \
-reducer NONE

$HADOOP_HOME/bin/hadoop fs -getmerge $outputbam $outputmergebam && $HADOOP_HOME/bin/hadoop fs -put $outputmergebam $outputmergebam |rm -rf $outputmergebam


$HADOOP_HOME/bin/hadoop jar $HADOOP_JAR \
-files hdfs://192.168.1.198:54310/user/hadoop/${Species},hdfs://192.168.1.198:54310/user/hadoop/$outputmergebam \
-input /user/hadoop/space.txt \
-output /user/hadoop/$dna_pair.sorted \
-mapper "$SAMTOOLS_HOME/samtools sort  $outputmergebam -o  $dna_pair.sorted " \
-reducer NONE



$SAMTOOLS_HOME/samtools rmdup $dna_pair.sorted.bam $outputrmbam

$HADOOP_HOME/bin/hadoop fs -put $dna_pair.sorted.bam $dna_pair.sorted.bam

$SAMTOOLS_HOME/samtools index $outputrmbam

$HADOOP_HOME/bin/hadoop fs -put $outputrmbam.bai $outputrmbam.bai


$HADOOP_HOME/bin/hadoop jar $HADOOP_JAR \
-files hdfs://192.168.1.198:54310/user/hadoop/${Species},hdfs://192.168.1.198:54310/user/hadoop/$outputrmbam,hdfs://192.168.1.198:54310/user/hadoop/$outputrmbam.bai,hdfs://192.168.1.198:54310/user/hadoop/${Species}.fai \
-input /user/hadoop/space.txt \
-output /user/hadoop/$outputrmbam.tobcf \
-mapper "$SAMTOOLS_HOME/samtools mpileup -C50 -q 30 -Bugf ${Species} $outputrmbam " \
-reducer NONE



$HADOOP_HOME/bin/hadoop jar $HADOOP_JAR \
-files hdfs://192.168.1.198:54310/user/hadoop/$outputrmbam.tobcf \
-input /user/hadoop/space.txt \
-output /user/hadoop/$outputrmbam.bcf \
-mapper "$BCFTOOLS/bcftools view -bcvg $outputrmbam.tobcf " \
-reducer NONE




$HADOOP_HOME/bin/hadoop jar $HADOOP_JAR \
-files hdfs://192.168.1.198:54310/user/hadoop/${Species},hdfs://192.168.1.198:54310/user/hadoop/$outputrmbam.bcf,hdfs://192.168.1.198:54310/user/hadoop/$outputrmbam.bai \
-input /user/hadoop/space.txt \
-output /user/hadoop/$outputrmbam.flt.vcf \
-mapper " $BCFTOOLS/bcftools view $outputrmbam.bcf | $BCFTOOLS/vcfutils.pl varFilter  " \
-reducer NONE