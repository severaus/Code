1.	导入数据. Importing data
2.	nohup gunzip -c file > newfile &   #解压gz文件，-c表示保留源文件
3.	nohup cat file1.txt file2.txt > file.txt &  #合并文件，每个样本的R1和R2分别合并到一块（每个样本的R1或者R2可能要测多次，将多次测得的R1合并到一块，多次测得的R2合并到一块），合并的时候合并顺序要一样。（将每个样本的R1和R2分别合并到一块后再合并为一个文件和直接将每个样本的所有测序序列直接合并为一个结果会一样吗？会影响后来的fastqc质控吗？ 答：R1和R2从来不需要合并到一块）
4.	FastQC #原始数据质控
nohup fastqc -f fastq -o 自建文件夹 file1 file2 file3 file4 &   #fastqc命令
5.	QC   ~24hff
Fastqc后由结果可以看出，后面有些reads结果不好，而且有些有很多N片段，所以需要去除，利用NGSQCToolkit_v2.3.3软件（http://blog.csdn.net/shmilyringpull/article/details/9225195）进行QC和trim。
首先利用IlluQC.pl 用于Illumina reads的QC。默认情况下去除掉含有primer/adaptor的reads和低质量的reads，并给出统计结果和6种图形结果。默认设置 (‘-s’ 参数) 碱基质量低于20的为低质量碱基；默认设置 ( ‘-l’ 参数)低质量碱基在reads中比例 >30% 的为低质量reads。本程序运行命令：
nohup perl /home/liul/software/NGSQCToolkit_v2.3.3/QC/IlluQC.pl -pe C1_R1.fastq C1_R2.fastq 2 A -p 2 -l 70 -s 20 -o C1_IlluQC &

nohup perl /home/liul/software/NGSQCToolkit_v2.3.3/QC/IlluQC.pl -pe F0C1_R1 F0C1_R2 2 A -p 4 -l 70 -s 20 -o F0C1_IlluQC &

再次，利用TrimmingReads.pl进行trim，（~2h）从3’端进行去掉碱基质量低于20的低质量碱基，去掉长度小于70的小片段。本程序运行命令：
nohup perl /home/liul/software/NGSQCToolkit_v2.3.3/Trimming/TrimmingReads.pl -i C1_R1.fastq_filtered -irev C1_R2.fastq_filtered -q 20 -n 70 -o C1_TrimmingReads &

nohup perl /home/liul/software/NGSQCToolkit_v2.3.3/Trimming/TrimmingReads.pl -i F0C1_R1_filtered -irev F0C1_R2_filtered -q 20 -n 70 -o F0C1_TrimmingReads &

QC建议用Trimmomatic进行，一步到位，省时间，这个软件（7806, bioinformatics）比NGS（1095，plos one）那个引用率高了6711次
nohup java -jar /mnt/disk0/leil/software/Trimmomatic-0-2.38/trimmomatic-0.38.jar PE -phred33 F0C1_1.fq.gz F0C1_2.fq.gz F0C1_R1.fq.gz F0C1_R1un.fq.gz F0C1_R2.fq.gz F0C1_R2un.fq.gz ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 &
(http://www.usadellab.org/cms/?page=trimmomatic)


6.	QC后用fastqc再跑下看看质控质量如何。

下面利用Hisat2和Stringtie等软件进行初步分析
7.	extract_exons.py Gallus.gtf> genome.exon  # 建立索引中的其中一步
8.	extract_splice_sites.py Gallus.gtf> genome.ss   # 建立索引中的其中一步
gtf和>之间不能加空格
#建立基因组+转录组索引：
bowtie2的索引只有基因组序列信息，tophat2比对时，转录组信息通过-G参数指定。HISAT2建立索引时，就应该把转录组信息加进去。
HISAT2提供两个Python脚本将GTF文件转换成hisat2-build能使用的文件。

9.	nohup hisat2-build –p 10 --ss genome.ss --exon genome.exon Gallus.fa gallus_tran &  #构建hisat2 index, 如果报错，可自己直接输入，因为可能是粘贴的时候有些中英字符不一样。
hisat2 -p 16 -x ./grch38_tran/genome_tran -1 SRR534293_1.fastq -2 SRR534293_2.fastq –S SRR534293.sam 
# -x 指定基因组索引   -1 指定第一个fastq文件   -2 指定第二个fastq文件  -S 指定输出的SAM文件
更多参数请查看HISAT2的操作手册：https://ccb.jhu.edu/software/hisat2/manual.shtml

比较全的hisat2-stringtie-ballgown的说明：http://blog.sciencenet.cn/home.php?mod=space&uid=1469385&do=blog&id=1022768
10.	nohup hisat2 –p 2 --dta –x gallus_tran -1 C1_R1 -2 C1_R2 –S C1.sam &  
#-p 线程； --x建立的参考基因组索引，gallus_tran是前缀，共有8个文件，-S 输出文件为sam文件。这里的代码是用于非特异性建库的，特异性建库的lncRNA需要用上--rna-strandness参数，如下所示:
 hisat2
•	非链特异性 默认，不需要设置
•	链特异性 如果使用dUTP: --rna-strandness RF， 添加了该参数后，每条 reads 将在 sam 文件中出现 XS 的tag，‘+’ ‘-’ 代表该 reads 所在的转录本与基因组序列的关系。
nohup hisat2 –p 2 --rna-strandness RF --dta –x gallus_tran -1 C1_R1 -2 C1_R2 –S C1.sam &

11.	nohup samtools sort -@ 8 -o C1.bam C1.sam &
#Sort and convert the SAM files to BAM 
12.	nohup stringtie -p 2 -G /mnt/disk0/leil/ref/Gallus.gtf –o C1.gtf –l C1 C1.bam & #Assemble transcripts for each sample （lnc分析再加上-f 0.01, -c 0.01. -f minimum isoform fraction (default: 0.1), -c minimum reads per bp coverage to consider for transcript assembly (default: 2.5)）
13.	nohup cuffcompare -i mergelist.txt -r /mnt/disk0/leil/ref/Gallus.gtf &
（https://blog.csdn.net/fengleqi/article/details/80623081）
#Merge transcripts from all samples.
  mergelist.txt内容
会输出4个总的文件，另外会在单个样品的gtf文件所在的位置每个样输出两个文件
Four total files will be output, and two files will be output for each sample at the location of the GTF file for a single sample.
换成cuffmerge更有优势（https://blog.csdn.net/shmilyringpull/article/details/8136837）Converting to cuffmerge is more advantageous
Cuffmerge
Create a file called assemblies.txt that lists the assembly file for each sample. The file should contain the following lines:
./C1_R1_clout/transcripts.gtf
./C2_R2_clout/transcripts.gtf
./C1_R2_clout/transcripts.gtf
./C2_R1_clout/transcripts.gtf-
./C1_R3_clout/transcripts.gtf
./C2_R3_clout/transcripts.gtf 
cuffmerge -o @@ -g genes.gtf -s genome.fa -p 8 assemblies.txt
nohup cuffmerge -o /diskd/liul/RNAseq_raw_data/cuffmerge -g /diskd/liul/refergenome/Gallus_gallus.Galgal4.83.gtf -s /diskd/liul/refergenome/Gallus.fa -p 15 /diskd/liul/RNAseq_raw_data/cuffmerge/assemblies.txt &


14.	perl summary_gtf.pl combined.gtf 200 2 length_exon.txt  
#利用Perl脚本summary_gtf.pl从combined.gtf 中筛选出长度大于等于200，外显子大于等于2的信息，输出文件是length_exon.txt
# Using the Perl script summary_gtf.pl, information with length greater than or equal to 200 and exon greater than or equal to 2 is screened out from combined.gtf. The output file is length_exon.txt.

15.	perl extract_gtf.pl length_exon.txt combined.gtf >length_exon.gtf
#利用length_exon.txt的信息从combined.gtf中提取信息，目的是将上一步筛选到的txt文件（length_exon.txt）文件转化为gtf文件（length_exon.gtf）
# Using the information of length_exon.txt to extract information from combined.gtf, the purpose is to transform the txt file filtered in the previous step (length_exon.txt) into GTF file (length_exon.gtf).

16.	grep –E ‘class_code “i”| class_code “u”| class_code “x”’ length_exon.gtf > iux.gtf     
 #从length_exon.gtf中筛选出有iux等的信息,输出文件是iux.gtf
# Screen out information such as iux from length_exon.gtf. The output file is iux.gtf.

grep ‘class_code “i”’ length_exon.gtf > i.gtf     
 #从length_exon.gtf中筛选出有i的信息,输出文件是i.gtf，筛选单个信息不需要加-E

Grep'class_code'i', length_exon.gtf > i.gtf

# Screen out the information with i from length_exon.gtf. The output file is i.gtf. Screening individual information does not need to add-E.


下面利用iux的gtf或者fa文件进行编码潜能预测
Next, use Iux GTF or FA files to predict coding potential
17.	编码潜能预测
17. Coding potential prediction
Transcripts with predicted protein-coding potential were removed (protein-coding potential criteria: CPC score > 0, PLEK score > 0, and CNCI score > 0). 
1）	CNCI （参考http://blog.biochen.com/archives/297及提到的软件链接）
Refer to http://blog.biochen.com/archives/297 and the software links mentioned

基本命令为：python CNCI_package/CNCI.py -f iux.fasta -o CNCI_out -m ve -p 10 #输出文件为CNCI.index和CNCI_out.log
The basic commands are: Python CNCI_package/CNCI.py-f iux.fasta-o CNCI_out-m ve-p 10  # output files are CNCI.index and CNCI_out.log.
参数说明：Description of parameters:
-f 输入fasta文件（可以使用-g参数输入GTF文件，但是同时需要使用
-d参数指定参考基因组的目录）
-o 输出结果目录
-m 指定模式，脊椎动物选择ve，植物选择pl
-p 指定CPU核数
- F Enter the FASTA file (you can enter the GTF file with the - G parameter, but you also need      to use the - D Parameter Designated Reference Genome Catalogue)
- O Output results directory
- M designated mode, vertebrate select ve, plant select pl
- P Specifies CPU Number

nohup python2 CNCI/filter_novel_lincRNA.py -i CNCI_out/CNCI.index -g /mnt/disk0/leil/lnc/8cuffcompare/filter/uix.gtf -s 0 -o out_dir &

这一步是紧接着上一步，可是输出Novel coding genes，Ambiguous genes，Novel lincRNA genes，Filter out noncoding genes等信息，但是只有lincRNA，不知道怎么整出其他种类的lncRNA，所以这一步暂时先不用，直接用上一步的CNCI.index进行后续分析
This step is followed by the next step, but the output of Novel coding genes, Ambiguous genes, Novel lincRNA genes, Filter out of non coding genes and other information, but only lincRNA, do not know how to integrate other types of lncRNA, so this step is not used for the time being, directly using the CNCI. index of the previous step for subsequent analysis.

2）	CPC（Coding Potential Calculator, 参考http://blog.biochen.com/archives/271 及提到的软件链接）
可以利用上面在线软件（http://cpc.cbi.pku.edu.cn/programs/run_cpc.jsp）直接跑，输入文件是iux.fa. 
2) CPC (Coding Potential Calculator, refer to http://blog.biochen.com/archives/271 and the software links mentioned)

The above online software (http://cpc.cbi.pku.edu.cn/programs/run_cpc.jsp) can be used to run directly. The input file is iux.fa.

3）	PLEK - predictor of lncRNAs and mRNAs based on k-mer scheme
python2 PLEK.1.2/PLEK.py -fasta /mnt/disk0/leil/lnc/8cuffcompare/filter/uix.gtf.fa -out PLEK_out -thread 12 -isoutmsg 1  #输出文件是PLEK_out，参数参考官网上的readme
3) PLEK-predictor of lncRNAs and mRNAs based on k-mer scheme

Python 2 PLEK.1.2/PLEK.py-fasta/mnt/disk0/leil/lnc/8cuffcompare/filter/uix.gtf.fa-out PLEK_out-thread 12-isoutmsg 1 output file is PLEK_out, parameters refer to readme on the official website.

4）	CPAT（参考http://blog.biochen.com/archives/219 及提到的软件链接）
4) CPAT (refer to http://blog.biochen.com/archives/219 and mentioned software links)
1) Download Gallus.cds.all.fa and Gallus.ncrna.fa from Ensemble
2) make_hexamer_tab.py -c Gallus.cds.all.fa -n Gallus.ncrna.fa > gallus_hexamer.tsv
    Output: gallus_hexamer.tsv
3) make_logitModel.py -x gallus_hexamer.tsv -c Gallus.cds.all.fa -n Gallus.ncrna.fa -o GAllus
    Output: GAllus.feature.xls, GAllus.logit.RData, GAllus.make_logitModel.r
4) cpat.py -g /mnt/disk0/leil/lnc/8cuffcompare/filter/uix.gtf.fa -d GAllus.logit.RData -x gallus_hexamer.tsv -o output
    Output files: output output.dat output.r
####coding probability的cutoff需要利用一个R程序进行判断，但是其中的train.dat不知道如何得到
#### The cutoff of coding probability needs to be judged by an R program, but train. dat does not know how to get it.

Human coding probability (CP) cutoff: 0.364 (CP >=0.364 indicates coding sequence, CP < 0.364 indicates noncoding sequence); Mouse cutoff: 0.44；Fly cutoff: 0.39；Zebrafish cutoff: 0.38
5）	筛选标准5) Screening criteria
##Transcripts with predicted protein-coding potential were removed (protein-coding potential criteria: CPC score > 0, PLEK score > 0, and CNCI score > 0). 


18.	
awk '{j=0;for(i=5;i<=NF;i++)if($i=="-")j++}{if(j<=7)print $1}' cuffcmp.tracking>atleast2.txt
#从cuffcmp.tracking里面筛选两个样均有的转录本
awk'{j=0; for (i=5; i<=NF; i++) if ($i=="-") j++}{if (j<=7) print $1}'cuffcmp. tracking > atleast2.txt

# Screening two transcripts from cuffcmp.tracking


19.	将18中的输出文件atleast2.txt和17中的几种软件预测的结果分别求交集，最后再将各交集进行求交集，得到几种编码潜能预测后的并在至少两个样本中有的转录本，输出文件保存为common.txt
#引入至少两个样本中有的转录本，是为了使找到的结果更准确，也缩小了数据量，利用后续分析
19. The output files atleast2.txt in 18 and 17 are intersected separately. Finally, the intersections are intersected to obtain transcripts predicted by several coding potentials and in at least two samples. The output files are saved as common txt.

# Transcripts from at least two samples were introduced to make the results more accurate and to reduce the amount of data. Follow-up analysis was used.

20.	perl extract_gtf.pl common.txt combined.gtf >common.gtf
#利用common.txt的信息从combined.gtf中提取信息，目的是将上一步筛选到的txt文件（common.txt）文件转化为gtf文件（common.gtf）进行后续差异表达的分析
awk -F "\"" '{print $4}' common.gtf|uniq|wc -l
#看gtf文件中有多少个转录本
wc -l common.gtf
#看gtf文件中有多少行

# Using common. TXT information to extract information from combined. gtf, the purpose is to transform the txt file filtered in the previous step into GTF file (common. gtf) for subsequent differential expression analysis.

Awk-F """'{print $4}'common.gtf | uniq | wc-l
# See how many transcripts are in the GTF file
Wc-l common.gtf
# See how many lines are in the GTF file


21.	nohup cuffdiff -o LC -b /mnt/disk0/leil/ref/Gallus.fa -p 15 -L L,C -u /mnt/disk0/leil/lnc/8cuffcompare/cuffcmp_output/common.gtf L1.bam,L2.bam,L3.bam C1.bam,C2.bam,C3.bam &
#计算差异表达的lncRNA，和普通转录组代码一样，只是gtf注释文件变为了筛选得到的有关lncRNA的gtf文件。注意：-L的每个组别标签用逗号隔开但不空格，每组的不同样本之间也是用逗号隔开但不空格，不同的分组之间用空格隔开,以上输出的结果三组在一起。两两之间做才能两两之间地输出结果

# Computing differentially expressed lncRNAs, like ordinary transcriptome codes, is only a GTF annotation file transformed into a screened GTF file about lncRNA. Note: Each group label of - L is separated by commas but not spaces, and different samples of each group are separated by commas but not spaces. Different groups are separated by spaces. The results of the above output are three groups together. Doing between two can output results between two.

22.	cis作用靶基因预测基本原理认为lncRNA的功能与其坐标临近的蛋白编码基因相关，于是将lncRNA临近位置的(上下游10k\100k)蛋白编码基因筛选出来作为其靶基因。后续再通过靶基因功能富集分析预测lncRNA的主要功能。

22. The basic principle of predicting target genes of CIS is that the function of lncRNA is related to the protein-coding genes near its coordinates. Therefore, the protein-coding genes near the location of lncRNA (upstream and downstream 10k 100k) are screened out as its target genes. The main functions of lncRNA were predicted by functional enrichment analysis of target genes.

将cuffdiff输出的差异表达lncRNA整理成如下格式：
The differentially expressed lncRNA of cuffdiff output was arranged in the following format:

 

(注意：只有转化为Windows Formmated Text才能用下面的指令成功转换 )
用sort -k1,1 -k2,2n LC.txt >LC.bed 将文件转化为bed文件。
用ensemble中的biomart将下载相应的注释信息
最后用bedtools进行寻找目标基因。
Note: Only when converted to Windows Formmated Text can the following instructions be used for successful conversion)
Use sort-k1, 1-k2, 2n LC.txt > LC.bed  to convert files into bed files.
The corresponding annotation information will be downloaded with biomaRt in ensemble
Finally, bedtools were used to search for target genes.

   bedtools window -a PL.bed -b Allgene.bed -w 100000 
   #-w表示寻找上下游位置，一般寻找10kb或者100kb。当时结果直接出现在指令编辑框内，复制粘贴出来的，无输出文件。
#- w means looking for upstream and downstream locations, usually 10 KB or 100 kb. At that time, the results appeared directly in the instruction edit box, copy and paste out, no output file.


23.	trans作用预测
差异表达的lncRNA和差异表达的mRNA进行相关性分析，找到trans预测的靶基因
23. Transaction prediction
The correlation between differentially expressed lncRNA and differentially expressed RNA was analyzed to find target genes for trans prediction.
setwd('/Users/leiliu/Desktop/')
lncRNA<-read.csv("PCDEL.csv",head=T)
mRNA<-read.csv("PCDEG.csv",head=T)
dim(lncRNA)
dim(mRNA)
head(lncRNA)
head(mRNA)
row.names(lncRNA) <- lncRNA[,1]
lncRNA <-lncRNA[, -1]
row.names(mRNA) <- mRNA[,1]
mRNA <- mRNA[, -1]
###### Cor ######
cor_matrix<-cor(t(lncRNA), t(mRNA),method = "spearman")
cor_matrix
install.packages("reshape")
library(reshape)
data <- melt(cor_matrix)
write.csv(data,"lncRNA_mRNA.correlation.spearman.csv")

###### P-value #####
install.packages("Hmisc")
library(Hmisc)
corpvalue <- rcorr(t(lncRNA), t(mRNA), type="spearman")
data <- melt(corpvalue$P)
write.csv(data,"lncRNA_mRNA.correlation.spearman.pvalue.csv")
由于相关系数和P值不能同时得到，利用两部，分别得到两个文件
Because correlation coefficient and P value can not be obtained simultaneously, two files can be obtained by using two parts.
利用linux处理cor和Pvalue的三个脚本，筛选出同时满足cor>0.9和P<0.05的条目
Using Linux to process three scripts of cor and Value, the entries satisfying both cor > 0.9 and P < 0.05 were screened out.

sed 's/\"//g' lncRNA_mRNA.correlation.spearman0.9.txt |awk '{print $2"|"$3"\t"$4}'|sort>cor.txt

sed 's/\"//g' lncRNA_mRNA.correlation.spearman.pvalue0.05.txt |awk '{if($2~/^XLOC/&& $3~/^ENS/)print $2"|"$3"\t"$4}'|sort>pvalue.txt

join cor.txt pvalue.txt >PL.txt




24.	利用预测的cis和trans靶基因的合集进行功能富集分析
24. Functional enrichment analysis using predicted aggregation of CIS and trans target genes
DAVID, KOBAS, IPA, STRING, Cytoscape
25.	sds
26.	fsds
27.	fsfs



 



28.	 利用awk进行初步筛选：保留FPKM>=0.5, length>200, exon number>=2的
28. Preliminary screening using awk: retain FPKM > = 0.5, length > 200, exon number > = 2
awk '{if($7>=0.5 && $10 > 1 && $11 >200) print $0}' cufcomp.OC_1yrF.stringtie.gtf.tmap > filter.OC_1yrF # $后面的数字代表相应选项的列数
Awk'{if ($7 >= 0.5 & & &$10 > 1 & &$11 > 200) print $0}'cufcomp. OC_1yrF. stringtie. gtf. TMAP > filter. OC_1yrF #$the following number represents the column number of the corresponding options

 awk '{if($6>=2 && $7 >=0.5 && $10 >200) print $0}' merged.stringtie_merged.gtf.tmap > filter.tmap1   
将下面命令加入环境变量
Awk'{if ($6 >= 2 & & & $7 >= 0.5 & & & $10 > 200) print $0}'merged. stringtie_merged. gtf. TMAP > filter. tmap1
Add the following command to the environment variable
export PATH=/path/to/install/hmmer3/bin:$PATH
export PERL5LIB=/path/to/pfam_scanDir:$PERL5LIB


 
extract_gtf.pl：
#! /usr/bin/perl 
#ZZ,YZ,YS,OZ,MS,DL,DB,CB
						
###################################################################

	open (T1,"$ARGV[0]");
	while($line=<T1>){
	chomp $line;
	@A =split (/\s+/ ,$line);
	$standrad{$A[0]}=1;
	}

	open (T2,"$ARGV[1]");
	while(<T2>){
	chomp;
	@A =split (/\"/ ,$_);
	if($standrad{$A[3]}==1){
	print "$_\n";
	}

	}	

summary_gtf.pl：
#!/usr/bin/perl -w
###########################################################################
# summary_gtf.pl --- Process gtf files
# <transcript_id> <gene_id> <chr> <start> <end> <len> <# of exons> <len_exon1;len_exon2;...>
###########################################################################
# Author: Fei Zhan <fei@fei-laptop>
# Created: 19 Dec 2012
# Version: 0.01
# Modified by Wentao at 22 Oct. 2018

use warnings;
use strict;

my $filein = "$ARGV[0]";
my $fileout = "$ARGV[3]";

my $min_length = "$ARGV[1]";
my $min_exons = "$ARGV[2]";

my $usage = "\nUsage:\n\t perl summary_gtf.pl -i [filein.gtf] -o [fileout.txt]\n\n";
#foreach my $i (0 ..scalar(@ARGV)-1) {
 # if($ARGV[$i] eq '-i') {
  #  $filein = $ARGV[++$i];
  #}elsif($ARGV[$i] eq '-o') {
   # $fileout = $ARGV[++$i];
  #}
#}
if(@ARGV ==0) {
    die $usage;
}

open(IN,'<',$filein) or die "Cannot open $filein.\n$usage";
open(OUT,'>',$fileout) or die "Cannot open $fileout.\n$usage";

my $headline = join("\t", ("TCONS", "XLOC", "chr", "strand", "start", "end", "num_exons", "length", "starts", "ends"));
print OUT $headline,"\n";

my %hash_gtf;
while(<IN>) {
  chomp;
  my @terms = split("\t",$_);
  my($chr,$start,$end,$strand,$info) = @terms[0,3,4,6,8];
  my $gene_id = $1 if $info =~ /gene_id "(.+?)";/;
  my $transcript_id = $1 if $info =~ /transcript_id "(.+?)";/;
  my $len_exon = $end-$start+1;
  if(exists $hash_gtf{$transcript_id}) {
    my $ref = $hash_gtf{$transcript_id};
    if($transcript_id ne $$ref[0] or $gene_id ne $$ref[1] or $chr ne $$ref[2] or $strand ne $$ref[3]) {
      print $_;
      die "Error!\n";
    }
    push(@{$$ref[4]}, $start);
    push(@{$$ref[5]}, $end);
    push(@{$$ref[6]}, $len_exon);
  }else{
    my @lens = ($len_exon);
    my @starts = ($start);
    my @ends = ($end);
    my $data = [$transcript_id,$gene_id,$chr,$strand,\@starts,\@ends,\@lens];
    $hash_gtf{$transcript_id} = $data;
  }
}

foreach my $key (keys %hash_gtf) {
  my $ref = $hash_gtf{$key};
  my @starts = @{$$ref[4]};
  my @ends = @{$$ref[5]};
  my $num_exons = scalar(@{$$ref[6]});
  my $len = 0;
  foreach (@{$$ref[6]}) {
    $len += scalar($_);
  }
  if($len >= $min_length && $num_exons >= $min_exons) {
    print OUT join("\t", (@{$hash_gtf{$key}}[0,1,2,3], $starts[0], $ends[-1], $num_exons, $len, join(",", @starts), join(",", @ends))),"\n";
  }
#  print OUT join("\t", (@{$hash_gtf{$key}}[0,1,2,3], $len, $num_exons)),"\n";
}


close IN;
close OUT;



__END__

=head1 NAME

process_gtf.pl - Describe the usage of script briefly

=head1 SYNOPSIS

process_gtf.pl [options] args

      -opt --long      Option description

=head1 DESCRIPTION

Stub documentation for process_gtf.pl, 

=head1 AUTHOR

Fei Zhan, E<lt>fei@fei-laptopE<gt>

=head1 COPYRIGHT AND LICENSE

Copyright (C) 2012 by Fei Zhan

This program is free software; you can redistribute it and/or modify
it under the same terms as Perl itself, either Perl version 5.8.2 or,
at your option, any later version of Perl 5 you may have available.

=head1 BUGS

None reported... yet.

=cut
 
如何抓取并查看gtf文件里的RNA类型
How to Grab and View RNA Types in GTF Files
leiliu@animalsci-138:/mnt/disk0/leil/ref$ awk -F ";" '{for(i=1;i<=NF;i++)if($i~/^gene_biotype/)print $i}' Gallus.gtf|head
leiliu@animalsci-138:/mnt/disk0/leil/ref$ awk -F ";" '{for(i=1;i<=NF;i++)if($i~/^ gene_biotype/)print $i}' Gallus.gtf|head
 gene_biotype "protein_coding"
 gene_biotype "protein_coding"
 gene_biotype "protein_coding"
 gene_biotype "protein_coding"
 gene_biotype "protein_coding"
 gene_biotype "protein_coding"
 gene_biotype "protein_coding"
 gene_biotype "protein_coding"
 gene_biotype "protein_coding"
 gene_biotype "protein_coding"
leiliu@animalsci-138:/mnt/disk0/leil/ref$ awk -F ";" '{for(i=1;i<=NF;i++)if($i~/^ gene_biotype/)print $i}' Gallus.gtf|wc -l
782749
leiliu@animalsci-138:/mnt/disk0/leil/ref$ wc -l Gallus.gtf
782754 Gallus.gtf
leiliu@animalsci-138:/mnt/disk0/leil/ref$ awk -F ";" '{for(i=1;i<=NF;i++)if($i~/^ gene_biotype/)print $i}' Gallus.gtf|sort|uniq -c
  26894  gene_biotype "lincRNA"
   3348  gene_biotype "miRNA"
    432  gene_biotype "misc_RNA"
      6  gene_biotype "Mt_rRNA"
     66  gene_biotype "Mt_tRNA"
 750090  gene_biotype "protein_coding"
    212  gene_biotype "pseudogene"
      6  gene_biotype "ribozyme"
    609  gene_biotype "rRNA"
     48  gene_biotype "scaRNA"
    699  gene_biotype "snoRNA"
    336  gene_biotype "snRNA"
      3  gene_biotype "sRNA"
leiliu@animalsci-138:/mnt/disk0/leil/ref$ awk '{if($3~/^transcript/)print}'|awk -F ";" '{for(i=1;i<=NF;i++)if($i~/^ gene_biotype/)print $i}' Gallus.gtf|sort|uniq -c
  26894  gene_biotype "lincRNA"
   3348  gene_biotype "miRNA"
    432  gene_biotype "misc_RNA"
      6  gene_biotype "Mt_rRNA"
     66  gene_biotype "Mt_tRNA"
 750090  gene_biotype "protein_coding"
    212  gene_biotype "pseudogene"
      6  gene_biotype "ribozyme"
    609  gene_biotype "rRNA"
     48  gene_biotype "scaRNA"
    699  gene_biotype "snoRNA"
    336  gene_biotype "snRNA"
      3  gene_biotype "sRNA"
^C
leiliu@animalsci-138:/mnt/disk0/leil/ref$ awk '{if($3~/^transcript/)print}' Gallus.gtf|awk -F ";" '{for(i=1;i<=NF;i++)if($i~/^ gene_biotype/)print $i}'|sort|uniq -c
   5972  gene_biotype "lincRNA"
   1116  gene_biotype "miRNA"
    144  gene_biotype "misc_RNA"
      2  gene_biotype "Mt_rRNA"
     22  gene_biotype "Mt_tRNA"
  30252  gene_biotype "protein_coding"
     43  gene_biotype "pseudogene"
      2  gene_biotype "ribozyme"
    203  gene_biotype "rRNA"
     16  gene_biotype "scaRNA"
    233  gene_biotype "snoRNA"
    112  gene_biotype "snRNA"
      1  gene_biotype "sRNA"
leiliu@animalsci-138:/mnt/disk0/leil/ref$ awk '{if($3~/^transcript/)print}' Gallus.gtf|awk -F ";" '{for(i=1;i<=NF;i++)if($i~/^ gene_b^Ctype/)print $i}'|sort|uniq -c
leiliu@animalsci-138:/mnt/disk0/leil/ref$ awk -F ";" '{for(i=1;i<=NF;i++)if($i~/^ Gallus.gtf^Cprint $i}' Gallus.gtf
leiliu@animalsci-138:/mnt/disk0/leil/ref$ grep -v 'protein_coding' protein_coding|head
grep: protein_coding: No such file or directory
leiliu@animalsci-138:/mnt/disk0/leil/ref$ grep -v 'protein_coding' Gallus.gtf|head
#!genome-build Gallus_gallus-5.0
#!genome-version Gallus_gallus-5.0
#!genome-date 2015-12
#!genome-build-accession NCBI:GCA_000002315.3
#!genebuild-last-updated 2016-12
1	ensembl	gene	292474	292583	.	+	.	gene_id "ENSGALG00000034440"; gene_version "1"; gene_name "gga-mir-6708"; gene_source "ensembl"; gene_biotype "miRNA";
1	ensembl	transcript	292474	292583	.	+	.	gene_id "ENSGALG00000034440"; gene_version "1"; transcript_id "ENSGALT00000076498"; transcript_version "1"; gene_name "gga-mir-6708"; gene_source "ensembl"; gene_biotype "miRNA"; transcript_name "gga-mir-6708-201"; transcript_source "ensembl"; transcript_biotype "miRNA";
1	ensembl	exon	292474	292583	.	+	.	gene_id "ENSGALG00000034440"; gene_version "1"; transcript_id "ENSGALT00000076498"; transcript_version "1"; exon_number "1"; gene_name "gga-mir-6708"; gene_source "ensembl"; gene_biotype "miRNA"; transcript_name "gga-mir-6708-201"; transcript_source "ensembl"; transcript_biotype "miRNA"; exon_id "ENSGALE00000372077"; exon_version "1";
1	ensembl	gene	317086	317182	.	+	.	gene_id "ENSGALG00000038707"; gene_version "1"; gene_name "gga-mir-1604"; gene_source "ensembl"; gene_biotype "miRNA";
1	ensembl	transcript	317086	317182	.	+	.	gene_id "ENSGALG00000038707"; gene_version "1"; transcript_id "ENSGALT00000061579"; transcript_version "1"; gene_name "gga-mir-1604"; gene_source "ensembl"; gene_biotype "miRNA"; transcript_name "gga-mir-1604-201"; transcript_source "ensembl"; transcript_biotype "miRNA";
leiliu@animalsci-138:/mnt/disk0/leil/ref$ grep -v 'protein_coding' Gallus.gtf>noncoding.gtf
leiliu@animalsci-138:/mnt/disk0/leil/ref$ 

备注：有个别code可能不对，可尝试着来
Note: Individual codes may be incorrect, but try it.

ls -ltrh Gallus.gtf 
#查看Gallus.gtf 文件大小

# View Gallus. GTF file size



 
背景：Background:
cis作用靶基因预测基本原理认为lncRNA的功能与其坐标临近的蛋白编码基因相关，于是将lncRNA临近位置的(上下游10k\100k)蛋白编码基因筛选出来作为其靶基因。后续再通过靶基因功能富集分析预测lncRNA的主要功能。
The basic principle of predicting target genes of CIS is that the function of lncRNA is related to the protein coding genes near its coordinates. Therefore, the protein coding genes near the location of lncRNA (up and down 10k 100k) are screened out as the target genes. The main functions of lncRNA were predicted by functional enrichment analysis of target genes.

求问：
我有了已有编码基因的染色体坐标信息，也有lncRNA的转录本坐标信息，要找出其对应的“上下游10k\100k”的蛋白编码基因应该怎么实现了？目前我想到的是读取每个lncRNA的上下游坐标信息，然后遍历基因信息输出，这样的话运算量有点大，循环太多遍了。感觉算法不太好。
有没有高手指点一二！？感激不尽！
Questions:
I have the chromosome coordinate information of the coding genes and the transcript coordinate information of lncRNA. How can I find out the corresponding protein coding genes of "up and down 10k\\\ 100k"? At present, what I think of is to read the coordinate information of each lncRNA upstream and downstream, and then traverse the output of gene information, which is a bit too much computation and too many cycles. The feeling algorithm is not very good.

Do you have a good finger to point out one or two!? Be deeply grateful!

答：Answer
这个最快捷的方法是用bedtools这个工具
将基因和lncRNA的位置储存在bed格式
bedtool里有个closedBed子程序可以找到最近的基因
使用windowBed可以找上下游特定距离的程序
希望对你有帮助

The quickest way to do this is to use bedtools

Store genes and lncRNA in bed format

There's a closedBed subroutine in bedtool to find the nearest gene

Programs that use Windows Bed to find specific distances upstream and downstream

I hope it will help you.

sort -k1,1 -k2,2n LC.txt >LC.bed

#将文件整理成染色体，起始位置，终止位置和基因名后，用sort排序，并转化为bed文件
# After sorting the files into chromosomes, starting positions, terminating positions and gene names, they are sorted and converted into bed files.

awk '{print $1"\t"$10"\t"$12"\t"$22"\t"$24}' common.gtf >extract.txt
#提取common.gtf中的第1，10，12，22，24列并用空格隔开（”\t”），输出文件为
# Extract columns 1, 10, 12, 22, 24 from common. GTF and separate them with spaces (" t"), and output the file as follows

extract.txt
wc –l AAA.txt     #查看AAA.txt文件有多少行
WC L AAA. TXT See how many lines are in the AAA. TXT file

awk < common.gtf '{print $1}' | sort | uniq -c
#查看common.gtf文件第一列内容的分类数及其对应出现的次数
# View the number of categories in the first column of the common. GTF file and the corresponding occurrences


awk '{ if ($4 >= 1 && $4 <= 10) print $1 }' sample.txt


rename "s/.mapped.ILLUMINA.bwa.CEU.low_coverage.20111****14.bam//" *
#名字太复杂的时候，进行批量重命名# Batch renaming when names are too complex
       

grep ‘lncRNA’ GRCg6a.gtf>lncRNA.gtf
#从GRCg6a.gtf中筛选出lncRNA并输出到lncRNA文件中
# Screening lncRNA from GRCg6a.gtf and exporting it to lncRNA file
