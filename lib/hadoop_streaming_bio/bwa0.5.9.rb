#!/usr/bin/ruby
#coding :utf-8
#require 'open3'

class  BWA0_5_9
 
  SPECIES = "dna_71_5000.fa"
  Data_Pair_1 = "testdata_1_1"
  Data_Pair_2 = "testdata_2_1"
  Data = "testdata"
  OutSam = "#{Data}.sam"
  Source1 = "http://192.168.1.198/hadoop/testdata/#{Data_Pair_1}.fq.gz"
  Source2 = "http://192.168.1.198/hadoop/testdata/#{Data_Pair_2}.fq.gz"
  Reference = "http://192.168.1.198/hadoop/#{SPECIES}"
  HADOOP_JAR= "/usr/local/hadoop/hadoop/contrib/streaming/hadoop-streaming-1.1.2.jar"
    
  def index
  
  
    getData = '~/Data/getData'
   
    #system "mkdir -p #{getData}|wget -p #{getData} #{Reference}"
   
    %x{$BAM_HOME/bwa index #{SPECIES}}

    %x{$HADOOP_HOME/bin/hadoop fs -put #{SPECIES}.amb  #{SPECIES}.amb}
    %x{$HADOOP_HOME/bin/hadoop fs -put #{SPECIES}.ann  #{SPECIES}.ann}
    %x{$HADOOP_HOME/bin/hadoop fs -put #{SPECIES}.bwt  #{SPECIES}.bwt}
    %x{$HADOOP_HOME/bin/hadoop fs -put #{SPECIES}.pac  #{SPECIES}.pac}
    %x{$HADOOP_HOME/bin/hadoop fs -put #{SPECIES}.rbwt #{SPECIES}.rbwt}
    %x{$HADOOP_HOME/bin/hadoop fs -put #{SPECIES}.rpac #{SPECIES}.rpac}
    %x{$HADOOP_HOME/bin/hadoop fs -put #{SPECIES}.rsa  #{SPECIES}.rsa}
    %x{$HADOOP_HOME/bin/hadoop fs -put #{SPECIES}.sa   #{SPECIES}.sa}
    
    system "rm -rf #{getData}/*"
  
  end
  
  def pre_process
    setData = '~/Data/setData/human_case'  
    #system "mkdir -p #{setData}|wget -p #{setData}/human_case   #{Source1}"
    #system "mkdir -p #{setData}|cd #{setData} | wget  #{Source1}"
    #system "mkdir -p #{setData}|cd #{setData} | wget  #{Source2}"
    system "gunzip #{setData}/*"
    %x{$HADOOP_HOME/bin/hadoop fs -put  #{setData}/#{Data_Pair_1}.fq  #{Data_Pair_1}.fq}
    %x{$HADOOP_HOME/bin/hadoop fs -put  #{setData}/#{Data_Pair_2}.fq  #{Data_Pair_2}.fq}
    system "touch spcace.txt|$HADOOP_HOME/bin/hadoop fs -put spcace.txt space.txt|rm -rf space.txt"
  end
  
  def aln_pair
    %x{$HADOOP_HOME/bin/hadoop jar #{HADOOP_JAR} \
    -files hdfs://192.168.1.198:54310/user/hadoop/#{Data_Pair_1}.fq,hdfs://192.168.1.198:54310/user/hadoop/#{SPECIES},hdfs://192.168.1.198:54310/user/hadoop/#{SPECIES}.bwt,hdfs://192.168.1.198:54310/user/hadoop/#{SPECIES}.rbwt  \
    -input /user/hadoop/space.txt \
    -output "#{Data_Pair_1}.sai" \
    -mapper "$BAM_HOME/bwa aln #{SPECIES} #{Data_Pair_1}.fq"  \
    -reducer NONE}
    
    %x{$HADOOP_HOME/bin/hadoop jar #{HADOOP_JAR} \
    -files hdfs://192.168.1.198:54310/user/hadoop/#{Data_Pair_2}.fq,hdfs://192.168.1.198:54310/user/hadoop/#{SPECIES},hdfs://192.168.1.198:54310/user/hadoop/#{SPECIES}.bwt,hdfs://192.168.1.198:54310/user/hadoop/#{SPECIES}.rbwt  \
    -input /user/hadoop/space.txt \
    -output "#{Data_Pair_2}.sai" \
    -mapper "$BAM_HOME/bwa aln #{SPECIES} #{Data_Pair_2}.fq"  \
    -reducer NONE}
  end
  
  def sampe
    %x{$HADOOP_HOME/bin/hadoop jar #{HADOOP_JAR} \
    -files hdfs://192.168.1.198:54310/user/hadoop/#{Data_Pair_1}.fq,hdfs://192.168.1.198:54310/user/hadoop/#{Data_Pair_2}.fq,hdfs://192.168.1.198:54310/user/hadoop/#{SPECIES},hdfs://192.168.1.198:54310/user/hadoop/#{Data_Pair_1}.sai,hdfs://192.168.1.198:54310/user/hadoop/#{Data_Pair_2}.sai,hdfs://192.168.1.198:54310/user/hadoop/#{SPECIES}.ann,hdfs://192.168.1.198:54310/user/hadoop/#{SPECIES}.amb,hdfs://192.168.1.198:54310/user/hadoop/#{SPECIES}.pac,hdfs://192.168.1.198:54310/user/hadoop/#{SPECIES}.bwt,hdfs://192.168.1.198:54310/user/hadoop/#{SPECIES}.rbwt,hdfs://192.168.1.198:54310/user/hadoop/#{SPECIES}.sa,hdfs://192.168.1.198:54310/user/hadoop/#{SPECIES}.rsa \
    -input  /user/hadoop/space.txt \
    -output "#{OutSam}"  \
    -mapper "$BAM_HOME/bwa sampe #{SPECIES} #{Data_Pair_1}.sai #{Data_Pair_2}.sai #{Data_Pair_1}.fq #{Data_Pair_2}.fq  " \
    -reducer NONE}
  end
  
  def aln_short
    %x{$HADOOP_HOME/bin/hadoop jar #{HADOOP_JAR} \
    -files hdfs://192.168.1.198:54310/user/hadoop/#{Data}.fq,hdfs://192.168.1.198:54310/user/hadoop/#{SPECIES},hdfs://192.168.1.198:54310/user/hadoop/#{SPECIES}.bwt,hdfs://192.168.1.198:54310/user/hadoop/#{SPECIES}.rbwt  \
    -input /user/hadoop/space.txt \
    -output "#{Data}.sai" \
    -mapper "$BAM_HOME/bwa aln #{SPECIES} #{Data}.fq"  \
    -reducer NONE}
  end

  def samse
    %x{$HADOOP_HOME/bin/hadoop jar #{HADOOP_JAR} \
    -files hdfs://192.168.1.198:54310/user/hadoop/#{Data}.fq,hdfs://192.168.1.198:54310/user/hadoop/#{SPECIES},hdfs://192.168.1.198:54310/user/hadoop/#{Data}.sai,hdfs://192.168.1.198:54310/user/hadoop/#{SPECIES}.ann,hdfs://192.168.1.198:54310/user/hadoop/#{SPECIES}.amb,hdfs://192.168.1.198:54310/user/hadoop/#{SPECIES}.pac,hdfs://192.168.1.198:54310/user/hadoop/#{SPECIES}.bwt,hdfs://192.168.1.198:54310/user/hadoop/#{SPECIES}.rbwt,hdfs://192.168.1.198:54310/user/hadoop/#{SPECIES}.sa,hdfs://192.168.1.198:54310/user/hadoop/#{SPECIES}.rsa \
    -input  /user/hadoop/space.txt \
    -output "#{OutSam}_test"  \
    -mapper "$BAM_HOME/bwa samse #{SPECIES} #{Data}.sai ##{Data}.fq  " \
    -reducer NONE}
  end 

  def bwasw
    %x{$HADOOP_HOME/bin/hadoop jar #{HADOOP_JAR} \
    -files hdfs://192.168.1.198:54310/user/hadoop/#{Data}.fq,hdfs://192.168.1.198:54310/user/hadoop/#{SPECIES},hdfs://192.168.1.198:54310/user/hadoop/#{Data}.sai,hdfs://192.168.1.198:54310/user/hadoop/#{SPECIES}.ann,hdfs://192.168.1.198:54310/user/hadoop/#{SPECIES}.amb,hdfs://192.168.1.198:54310/user/hadoop/#{SPECIES}.pac,hdfs://192.168.1.198:54310/user/hadoop/#{SPECIES}.bwt,hdfs://192.168.1.198:54310/user/hadoop/#{SPECIES}.rbwt,hdfs://192.168.1.198:54310/user/hadoop/#{SPECIES}.sa,hdfs://192.168.1.198:54310/user/hadoop/#{SPECIES}.rsa \
    -input  /user/hadoop/space.txt \
    -output "#{OutSam}_test_2"  \
    -mapper "$BAM_HOME/bwa bwasw #{SPECIES} #{Data}.sai ##{Data}.fq  " \
    -reducer NONE}
  end

end
