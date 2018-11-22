# 根据文献 Label-free Prostate Cancer Detection by Characterization of Extracellular Vesicles using Raman Spectroscope 
来处理raman光谱数据，数据下机后用采用自带的软件进行去背景，归一化等，得到各样本波谱~强度表

# 寻找波长范围
文献中PCA选取的范围为 ![wave](https://github.com/xiao7462/Bioinformatics/blob/master/pic/wave.png) ，但我们的实验是不知道的，所以只有根据文献的方法遍历所有的可能性，找到一个最大的PCA区分图

# [代码记录](https://github.com/xiao7462/Bioinformatics/blob/master/raman_spectrum/pca2.py)
  笔记本带不动，在服务器上运行，得到170M遍历的图像  ![du-lh](https://github.com/xiao7462/Bioinformatics/blob/master/pic/du%20-lh.png)
