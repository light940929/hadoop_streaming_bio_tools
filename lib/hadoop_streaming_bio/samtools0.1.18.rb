#!/usr/bin/ruby
#coding :utf-8
require 'open3'

class Samtools0_1_18
  
  SPECIES = "dna_71_5000.fa"
  Data_Pair_1 = "testdata_1_1"
  Data_Pair_2 = "testdata_2_1"
  Data = "testdata"
  OutSam = "#{Data}.sam"
  OutBam = "#{Data}.bam"
  HADOOP_JAR= "/usr/local/hadoop/hadoop/contrib/streaming/hadoop-streaming-1.1.2.jar"
    
  %x{$HADOOP_HOME/bin/hadoop fs -getmerge #{OutSam} #{OutSam}_m && $HADOOP_HOME/bin/hadoop fs -put #{OutSam}_m #{OutSam}_m | rm -rf #{OutSam}_m}
  
  
  def faidx
    system "$SAMTOOLS_HOME/samtools faidx  #{SPECIES}"
    %x{$HADOOP_HOME/bin/hadoop fs -put #{SPECIES}.fai #{SPECIES}.fai | rm -rf #{SPECIES}.fai}
  end
   
  def toSam
    %x{$HADOOP_HOME/bin/hadoop jar #{HADOOP_JAR} \
    -files hdfs://192.168.1.198:54310/user/hadoop/#{SPECIES},hdfs://192.168.1.198:54310/user/hadoop/#{SPECIES}.fai,hdfs://192.168.1.198:54310/user/hadoop/#{OutSam}_m \
    -input /user/hadoop/space.txt \
    -output "/user/hadoop/#{OutBam}" \
    -mapper "$SAMTOOLS_HOME/samtools view -bt #{SPECIES}.fai #{OutSam}_m  \
    -reducer NONE}
  end
  
  def merge
    %x{$HADOOP_HOME/bin/hadoop fs -getmerge #{OutBam} #{OutBam}_m && $HADOOP_HOME/bin/hadoop fs -put #{OutBam}_m #{OutBam}_m|rm -rf #{OutBam}_m}
  end
  
  def sort  
    %x{$HADOOP_HOME/bin/hadoop jar #{HADOOP_JAR} \
    -files hdfs://192.168.1.198:54310/user/hadoop/#{SPECIES},hdfs://192.168.1.198:54310/user/hadoop/#{OutBam}_m \
    -input /user/hadoop/space.txt \
    -output "/user/hadoop/#{OutBam}_m.sorted" \
    -mapper "$SAMTOOLS_HOME/samtools sort  #{OutBam}_m -o  #{OutBam}_m.sorted " \
    -reducer NONE}
  end
  
  def rmdup 
    %x{$SAMTOOLS_HOME/samtools rmdup #{OutBam}_m.sorted.bam #{OutBam}_m_rmdup}
    %x{$HADOOP_HOME/bin/hadoop fs -put #{OutBam}_m_rmdup #{OutBam}_m_rmdup}
  end
  
  def index
    %x{$SAMTOOLS_HOME/samtools index #{OutBam}_m_rmdup}
    %x{$HADOOP_HOME/bin/hadoop fs -put #{OutBam}_m_rmdup.bai #{OutBam}_m_rmdup.bai}
  end
  
  def mpileup
  
    %x{$HADOOP_HOME/bin/hadoop jar #{HADOOP_JAR} \
    -files hdfs://192.168.1.198:54310/user/hadoop/#{SPECIES},hdfs://192.168.1.198:54310/user/hadoop/#{OutBam}_m_rmdup,hdfs://192.168.1.198:54310/user/hadoop/#{OutBam}_m_rmdup.bai,hdfs://192.168.1.198:54310/user/hadoop/#{SPECIES}.fai \
    -input /user/hadoop/space.txt \
    -output "/user/hadoop/#{OutBam}_m_rmdup}.tobcf"\
    -mapper "$SAMTOOLS_HOME/samtools mpileup -C50 -q 30 -Bugf #{SPECIES} #{OutBam}_m_rmdup" \
    -reducer NONE}
  
  end
  
  def bcfg
  
    %x{$HADOOP_HOME/bin/hadoop jar #{HADOOP_JAR} \
    -files hdfs://192.168.1.198:54310/user/hadoop/#{OutBam}_m_rmdup}.tobcf \
    -input /user/hadoop/space.txt \
    -output "/user/hadoop/#{OutBam}_m_rmdup.bcf" \
    -mapper "$BCFTOOLS/bcftools view -bcvg #{OutBam}_m_rmdup}.tobcf " \
    -reducer NONE}
  
  end
  
  def vcf
   
    %x{$HADOOP_HOME/bin/hadoop jar #{HADOOP_JAR} \
    -files hdfs://192.168.1.198:54310/user/hadoop/#{SPECIES},hdfs://192.168.1.198:54310/user/hadoop/#{OutBam}_m_rmdup.bcf,hdfs://192.168.1.198:54310/user/hadoop/#{OutBam}_m_rmdup.bai \
    -input /user/hadoop/space.txt \
    -output "/user/hadoop/#{OutBam}_m_rmdup.flt.vcf" \
    -mapper " $BCFTOOLS/bcftools view #{OutBam}_m_rmdup.bcf | $BCFTOOLS/vcfutils.pl varFilter  " \
    -reducer NONE}
  
  end
  
  
end