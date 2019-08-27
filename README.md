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
-[绘图](#绘图)

## 转录组
 * sra转fastq格式
     ```
     for i in SRR3159774 SRR3159776 SRR3159778 SRR3159775 SRR3159777 SRR3159779
     do 
      /data1/tangx/software/sratoolkit.2.9.2-ubuntu64/bin/fasterq-dump.2.9.2 --split-3 -e 20 -p ${i}
     done  #转换为fastq格式
     ```

## 无参转录组
  * 差异表达分析
  * [单独查看表达基因](https://github.com/xiao7462/Bioinformatics/blob/master/RNA-seq/%E5%AF%8C%E9%9B%86%E5%85%A8%E9%83%A8%E8%A1%A8%E8%BE%BE%E5%9F%BA%E5%9B%A0.md)

## 进化树

## 生信软件安装

## 拉曼光谱数据处理]

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
