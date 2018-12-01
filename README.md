# Bioinformatics
**This is a bioinfomation study note repository**
# RNA-seq
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
# Evolution tree

## 转录组进化树

# SNP

# Software

# 拉曼光谱PCA

# 代谢组处理
 * [数据处理](https://github.com/xiao7462/Bioinformatics/blob/master/metabolome/metabolome.ipynb) 待完成
