**Bioinformatics**     
**This is a bioinfomation study note repository**           
## 目录
-[转录组](#转录组)    
-[无参转录组](#无参转录组)         
-[进化树](#进化树)   
-[生信软件安装](#生信软件安装)   
-[拉曼光谱数据处理](#拉曼光谱数据处理)     
-[代谢组处理](#代谢组处理)    
-[PCR数据处理](#PCR数据处理)      
-[绘图](#绘图)<br>	
-[小脚本](#小脚本)
## 转录组
 * NCBI 数据下载
 ```bash
 for i in ` seq 56 62`;
do
    wget ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByStudy/sra/SRP/SRP075/SRP075747/SRR35899${i}/SRR35899${i}.sra  #也可以使用axel命令代替wget，axel比wget快一些
done
#56-58的舍弃了，只用 59-62
 ```
 
 * 通过sra下载ncbi数据
 ![sra.png](https://i.loli.net/2019/08/27/YBtNJEMLuDWeiRS.png)
 
 * 比对软件[总结](https://mp.weixin.qq.com/s/YI8QzAaAEWubCe1JxXEL1w?)
 

 * 质控+查看质控
 ```bash
 for i in `seq 56 62`
do 
    fastq-dump --gzip --split-3 -O ./fastq/ -A SRR35899${i}.sra
done  #转换为fastq格式

fastqc *   #质量检测

multiqc .  #一次性查看qc
 ```
 * 序列比对  **hisat2**
 ```bash
 for i in `seq 59 62`;
do
hisat2 -t -p 20 -x ~/reference/genome/mm10/mm10/genome\
 -1 ~/new_rnaseq/fastq/SRR35899${i}.sra_1.fastq.gz\
 -2 ~/new_rnaseq/fastq/SRR35899${i}.sra_2.fastq.gz\
 -S ~/new_rnaseq/bam/SRR35899${i}.sam >${i}.hisat2.log &&
samtools view -S SRR35899${i}.sam -b > SRR35899${i}.bam &&  #转换成bam文件
samtools sort SRR35899${i}.bam -o SRR35899${i}_sorted.bam &&  #采用默认方式排序
samtools index SRR35899${i}_sorted.bam  #生成索引
done
 ```
 
 * sra转fastq格式
     ```
     for i in SRR3159774 SRR3159776 SRR3159778 SRR3159775 SRR3159777 SRR3159779
     do 
      /data1/tangx/software/sratoolkit.2.9.2-ubuntu64/bin/fasterq-dump.2.9.2 --split-3 -e 20 -p ${i}
     done  #转换为fastq格式
     ```
* 下载imicrobe上所有的数据
```bash
for i in `seq 1662 2522`;
do
	wget ftp://ftp.imicrobe.us/projects/104/samples/${i}/*.nt.fa.gz
done

```

* sam转换为bam文件
```bash
sam--- bam  
samtools view -S siSUZ12_1.sam -b > siSUZ12_1.bam   
samtools sort siSUZ12_1.bam siSUZ12_1.sorted 
```

* read计数
代码
```python
！/usr/bin/python3
def combine_count(count_list):   ##定义一个合并的函数
    mydict={}                    #创建字典
    for file in count_list:       #读取count_list里面的每一个文件
        for line in open(file,'r'):      #读取每一个文件里面的每一行
            #print(line)
            if line.startswith('E'):      #判断是不是每行的开头都是E
                key,value = line.strip().split('\t')  #读取文件内容到新的字典里； \t 是制表符；  split('\t)是去除制表符；strip()是去除字符串头尾的空格以及 \n, \t之类的
                if key in mydict:
                    mydict[key]=mydict[key] + '\t' + value      # 追加值    
                else:ne
                    mydict[key]=value              # 更新key的内容
    sorted_mydict = sorted(mydict)      #排序
    out = open('count_out.txt','w')     #open创建一个新的空文件
    
    for k in sorted_mydict:              
        #print(k,mydict[k]）
        #break
        out.write(k+'\t'+mydict[k]+'\n')    # 写入已经排好序的内容
        
count_list=['SRR3589959.count', 'SRR3589960.count', 'SRR3589961.count'，'SRR3589962.count']

combine_count(count_list)

```

* Deseq2 差异表达分析
```r
library(DESeq2)
                ##数据预处理
                database <- read.table(file = "mouse_all_count.txt", sep = "\t", header = T, row.names = 1)
                database <- round(as.matrix(database))
                 
                ##设置分组信息并构建dds对象
                condition <- factor(c(rep("control",2),rep("Akap95",2)), levels = c("control", "Akap95"))
                coldata <- data.frame(row.names = colnames(database), condition)
                dds <- DESeqDataSetFromMatrix(countData=database, colData=coldata, design=~condition)
                dds <- dds[ rowSums(counts(dds)) > 1, ]
                 
                ##使用DESeq函数估计离散度，然后差异分析获得res对象
                dds <- DESeq(dds)
                res <- results(dds)
                 
                #最后设定阈值，筛选差异基因，导出数据(全部数据。包括标准化后的count数)
                res <- res[order(res$padj),]
                diff_gene <- subset(res, padj < 0.05 & (log2FoldChange > 1 | log2FoldChange < -1))
                diff_gene <- row.names(diff_gene)
                resdata <- merge(as.data.frame(res), as.data.frame(counts(dds, normalized=TRUE)),by="row.names",sort=FALSE)
                write.csv(resdata,file = "control_vs_Akap95.csv",row.names = F)

```

* Go 富集分析

```r
library(clusterProfiler)
                library(org.Mm.eg.db)
                 
                ##去除ID的版本号
                diff_gene_ENSEMBL <- unlist(lapply(diff_gene, function(x){strsplit(x, "\\.")[[1]][1]}))
                ##GOid mapping + GO富集
                ego <- enrichGO(gene = diff_gene_ENSEMBL, OrgDb = org.Mm.eg.db,
                                keyType = "ENSEMBL", ont = "BP", pAdjustMethod = "BH",
                                pvalueCutoff = 0.01, qvalueCutoff = 0.05)
                ##查看富集结果数据
                enrich_go <- as.data.frame(ego)
                ##作图
                barplot(ego, showCategory=10)
                dotplot(ego)
                enrichMap(ego)
                plotGOgraph(ego)

```
![pathway.png](https://i.loli.net/2019/08/27/u3NBe9hUXnSOFi7.png)


## 无参转录组
  * 差异表达分析
  * [单独查看表达基因](https://github.com/xiao7462/Bioinformatics/blob/master/RNA-seq/%E5%AF%8C%E9%9B%86%E5%85%A8%E9%83%A8%E8%A1%A8%E8%BE%BE%E5%9F%BA%E5%9B%A0.md)

## 进化树
### 进化树简单入门
[BLAST+ mafft +RAxml](http://www.biotrainee.com/thread-1521-1-1.html)

## 生信软件安装

### trimmomatics
实际运行代码
```bash
nohup java -jar /data1/tangx/software/Trimmomatic/Trimmomatic-0.36/trimmomatic-0.36.jar PE -phred33 ~/cnv/trim/A2-B2-3/A2-B2-3_1P.fq.gz ~/cnv/trim/A2-B2-3/A2-B2-3_2P.fq.gz -baseout ~/cnv/trim/A2-B2-3/A2-B2-3.fq.gz ILLUMINACLIP:/data1/tangx/software/Trimmomatic/Trimmomatic-0.36/adapters/TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 HEADCROP:8 MINLEN:36 & 
```
ILLUMINACLIP： 一定要注意下接头的路径更改
注意文件命名格式，带有不合规的符号要在路径上加 ""

批处理
```bash
for i in A01 A02 A03 A121 A122 A123 A41 A42 A43 A81 A82 A83;
do
java -jar /data1/tangx/software/Trimmomatic/Trimmomatic-0.36/trimmomatic-0.36.jar PE -phred33 ~/rna-seq/A/${i}.R1.fastq.gz ~/rna-seq/A/${i}.R2.fastq.gz -baseout ~/rna-seq/A/trim/${i}.fq.gz ILLUMINACLIP:/data1/tangx/software/Trimmomatic/Trimmomatic-0.36/adapters/TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 HEADCROP:8 MINLEN:36
done

##感觉还可以加个 && rm ~/rna-seq/A/${i}.R1.fastq.gz ~/rna-seq/A/${i}.R2.fastq.gz

## 继续用fastqc 连续查看  
for i in ~/rna-seq/A/trim/*P.fq.gz ;
do
    fastqc ${i} ;
done
```
[参数设置](https://www.jianshu.com/p/a8935adebaae)         

### Trinity
* 安装注意
> 直接用移动硬盘里面的安装包来安装trinity
> 一定先先安装bowtie2并且添加到~/.bashrc

* 服务器实际运行代码
```bash
nohup /data1/tangx/software/trinityrnaseq-Trinity-v2.4.0/Trinity --no_version_check --seqType fq --left ~/rna-seq/E-B/trim/E01/E01.fa_1P.gz --right ~/rna-seq/E-B/trim/E01/E01.fa_1P.gz  --CPU 30 --max_memory 100G --output E01_trinity &1>E01.log 2>E01.err &
```
* 注意点
下载ncbi上的数据注意点,trintiy需要特定的fastq-dump转换格式
```bash
fastq-dump --defline-seq '@$sn[_$rn]/$ri' --split-files file.sra --split-3 -O|--outdir  # 输出的文件为trinity

```
[参数设置](https://github.com/trinityrnaseq/trinityrnaseq/wiki/Running-Trinity)
一些教程
> [官网](https://github.com/trinityrnaseq/trinityrnaseq/wiki) [Trinity的安装与使用](http://blog.sina.com.cn/s/blog_90e51ee00102uwr0.html) [Trinity的安装与使用-生信技能树anlan](http://www.biotrainee.com/thread-2109-1-1.html) [Trinity进行转录组分析的一条龙服务](http://yangl.net/2016/12/16/trinity/)              


### samtools 
* 安装 
直接` sudo apt-get install samtools` 
[参数](https://blog.csdn.net/luobailian/article/details/50316627)

运行代码`samtools sort -@ 5 -n -m 20G ~/cnv/A2-B4-3/A2-B4-3.bam -o A2-B4-3.sorted-n.bam  # 注意 线程数和 内存数相乘不能超过服务器的配置
`


### bwa

安装 
```bash
git clone https://github.com/lh3/bwa.git
cd bwa; make
./bwa index ref.fa
./bwa mem ref.fa read-se.fq.gz | gzip -3 > aln-se.sam.gz
./bwa mem ref.fa read1.fq read2.fq | gzip -3 > aln-pe.sam.gz

```


### fastqc               
* 安装
FastQC官网

来自 <http://www.bioinformatics.babraham.ac.uk/projects/fastqc/> 

官网直接下载压缩文件安装，添加到~/.bashrc。   注意不要apt-get安装，会出bug

* 使用
fastqc 1.fa 2.fa

### cnvnator             
* 安装 
```bash
①CNVnator_v0.3.3.zip

来自 <https://github.com/abyzovlab/CNVnator/releases> 
下载

② wget https://root.cern.ch/download/root_v6.13.02.Linux-ubuntu16-x86_64-gcc5.4.tar.gz  #一定要下载Binary distributions  来下载

③tar zxvf root压缩包  
   cd root 
export ROOTSYS=/data1/tangx/software/root
export LD_LIBRARY_PATH=/data1/tangx/software/root/lib/root:${LD_LIBRARY_PATH}


④cd  CNVnator/src/samtools   
	Make 
	
Cd ..
make
```
* cern root 安装 
```bash
①下载安装包并解压 ， cd 进去 
②创建root6-build  ，cd 进去
③cmake ..
④cmake .   && make -j4>make.out 2>&1


安装cnvnator

$ cd src/samtools
$ make
Even if compilation is not completed, but the file libbam.a has been created, you
can continue.
$ cd ../
$ make 

```
参考教程
最近正在做这个，用一组有ground truth的数据尝试了3个软件：
CNVnator GitHub FP很高，难以做下游过滤
BIC-seq2 还在run...
CNVkit GitHub 比较成功，通过适当的参数选择可以得出sensitivity和specificity都是100%得结果

 https://www.zhihu.com/question/51602122/answer/297048446 

CNVnator Based Analysis

来自 <https://wlz0726.github.io/2017/04/27/CNVnator-Based-Analysis/> 

CNVnator Install & Usage
 

来自 <http://blog.sina.com.cn/s/blog_8fef6fbb0101nj2v.html> 

linux下安装ROOT过程

来自 <https://blog.csdn.net/i826056899/article/details/39596967> 

### RAxml
* 安装
![RAxml.png](https://i.loli.net/2019/08/27/9gan64iedsfYyxD.png)

* 参数
![parameter.png](https://i.loli.net/2019/08/27/lA1k6bhCHZR53Jq.png)

* 实际代码  -计算自展值(bootstrp)
`raxmlHPC-PTHREADS-AVX -T 4 -f a -p 12345 -x 12345 -# 100  -m PROTGAMMAWAG -sap2_align_cured.phy -n T9` 


### hisat2 
* 运行代码
```bash
extract_exons.py Homo_sapiens.GRCh38.83.chr.gtf > genome.exon
extract_splice_sites.py Homo_sapiens.GRCh38.83.chr.gtf > genome.ss
extract_snps.py snp142Common.txt > genome.snp
hisat2-build -p 4 genome.fa --snp genome.snp --ss genome.ss --exon genome.exon genome_snp_tran

```
### orthofinder
* 实际运行
`orthofinder -f ~/polytrans/raw_data -t 40`
`orthofinder -b ~/polytrans/Results_May20/WorkingDirectory -f ~/polytrans/E-B -t 40`



## 拉曼光谱数据处理]
```python
import matplotlib
matplotlib.use('Agg')
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os
from sklearn.decomposition import PCA
import sys

df = pd.read_csv('/data1/tangx/raman/0323/0323.csv',sep = ',')
df.set_index('Wave',inplace=True)
df = df.T
P_sample =df[df['target'] =='N' ]
Z_sample = df[df['target'] =='Z' ]
df = P_sample.append(Z_sample)
y = df['target'].values
X = df.iloc[:,1:]

wave_lengths = range(4,21)
for wave_length in wave_lengths:
    for i in range(500,1100): 
        wave_choose = X.iloc[:,i :i +wave_length]
        left = wave_choose.columns[0]
        right = wave_choose.columns[-1]
        # PCA
        pca = PCA(n_components = 3)
        principalComponents = pca.fit_transform(wave_choose)
        """principalComponents
        array([[ 0.11424775, -0.02477196,  0.00409868],
       [ 0.01258328,  0.05578577, -0.00444276],
       [-0.19815818,  0.03132211, -0.01661246],
       [ 0.05201658, -0.05212701, -0.00767984],
       [-0.08611526,  0.02718982,  0.02507499],
       [-0.01478776, -0.02974936, -0.01112425],
       [ 0.2506098 ,  0.04360269, -0.00337041],
       [-0.0476935 , -0.00448015, -0.00323041],
       [-0.01101285, -0.02548083,  0.01021661],
       [-0.07168986, -0.02129109,  0.00706986]])
        """
        principalDf = pd.DataFrame(data = principalComponents
             , columns = ['principal component 1', 'principal component 2','principal component 3'])
        principalDf['target'] = df['target'].values
        finalDf = principalDf
        
        # 绘图
        for pc_number in ([0,1],[1,2],[0,2]):
            plt.figure(figsize= (15,10))
            plt.xlabel('Principal Component {} ({:.2%})'.format(pc_number[0]+1,pca.explained_variance_ratio_[pc_number[0]]), fontsize = 15)
            plt.ylabel('Principal Component {} ({:.2%})'.format(pc_number[1]+1,pca.explained_variance_ratio_[pc_number[1]]), fontsize = 15)
            plt.title('{} to {} of raman wave'.format(left,right), fontsize = 20)
            targets = ['N','Z']
            colors = ['navy', 'red']
            for target, color in zip(targets,colors):
                indexs = finalDf[finalDf.target==target].index.tolist()
                plt.scatter(finalDf.iloc[indexs, pc_number[0]]
                           , finalDf.iloc[indexs, pc_number[1]]
                           , c = color,alpha=.8
                           , s = 50)
            plt.legend(targets)
            plt.grid()
            plt.savefig("./picture/{} to {} with {} to{}.png".format(right,left,pc_number[0]+1,pc_number[1]+1))
            #plt.show()
            plt.close()
```


## 代谢组处理
 * [数据处理](https://github.com/xiao7462/Bioinformatics/blob/master/metabolome/metabolome.ipynb) 待完成

## PCR数据处理
[q-PCR数据处理](https://github.com/xiao7462/Bioinformatics/blob/master/q-PCR.md)

## 绘图

### 热图绘制
 * 热图绘制 R     
 ![原始数据](https://github.com/xiao7462/Bioinformatics/blob/master/pic/heatmap_rawdata.png)
 ``` R
 rm(list=ls())
library(pheatmap) 
X = read.csv("C:/Users/tangxing/Desktop/111.csv",header = T)
rownames(X) <- X[,1] #将行名设置为列表的行名
X <- X[,-1] # 去掉第一行
X<-log(X,10) #  将数值转换
X[sapply(X,is.infinite)]<- -3 # 替换负无穷为-3
pheatmap(X,cellwidth = 30, cellheight = 9,border_color="black") # 绘制热图
 ```
 ![热图](https://github.com/xiao7462/Bioinformatics/blob/master/pic/heatmap.png)

* ggplot2绘制富集图
```r
library(ggplot2)
 
# 读取数据
pathway = read.table("./qwe.txt",header=T,sep="\t")
 
# 开始画图
p = ggplot(pathway,aes(richFactor,Pathway))
p + geom_point()
 
# 改变点的大小
p + geom_point(aes(size=Input-number))
 
# 四维数据的展示
pbubble = p + geom_point(aes(size=Input-number,color=-1*log10(Qvalue)))
 
# 自定义渐变颜色
pbubble + scale_colour_gradient(low="green",high="red")
 
# 绘制pathway富集散点图
pr = pbubble + scale_colour_gradient(low="green",high="red") + labs(color=expression(-log[10](Qvalue)),size="Gene number",x="Rich factor",y="Pathway name",title="Top20 of pathway enrichment")
# 改变图片的样式（主题）
pr + theme_bw()
## 保存图片
ggsave("out.pdf")   # 保存为pdf格式，支持 pdf，png，svg多重格式
ggsave("out.png")  # 保存为png格式
ggsave("out2.png",width=4,height=4)   # 设定图片大小
```

* 散点图
```r
library(ggplot2)
aa<-read.table("./444.txt",header=T,sep='\t')
ggplot(data=aa,aes(x=factor(batch),y=value,size=value,colour=factor(serial)))+
geom_point()+
#geom_smooth(method=lm,se=FALSE) +
#geom_boxplot(mapping = aes(group =factor(batch)), fill = 'steelblue')+
theme_gray()+
labs(x="batch",y="value")

```

## 小脚本
* 将基因按照随机的长度进行分隔 
```python
f = open('yelvti.txt')
lines = f.read()
lines = lines.replace('\n',''); lines
fw = open('random700-1300BP.fasta','w')
prev = 0
while True:
    n = random.randint(700,1300)
    # n = 1000时是固定为1000BP
    #splitted.append(lines[0][prev:prev+n])
    fw.write('>{} to {}'.format(prev, prev+n) + '\n')
    fw.write(lines[0][prev:prev+n] + '\n')
    prev += n 
    if prev > len(lines[0] ) - 1:
        break
```

[衣藻参考基因组](ftp://ftp.ensemblgenomes.org/pub/plants/release-44/fasta/chlamydomonas_reinhardtii/dna/)

