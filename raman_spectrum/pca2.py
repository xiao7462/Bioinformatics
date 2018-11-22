import matplotlib
matplotlib.use('Agg') # linux无GUI接口时
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os
from sklearn.decomposition import PCA
import sys
#np.set_printoptions(threshold='nan')

#数据预处理：加载数据，添加label
df = pd.read_csv(sys.argv[1],sep = '\t')
df.set_index('Wave',inplace=True)
X = df.T
y = np.zeros(200)
y[:49]=0
y[49:99]=1
y[99:150]=2
y[150:]=3

#去除异常值
qwe = X
qwe['y'] = y
qwe.head()
p1 = qwe[y==0]
N1  = qwe[y==1]
p0 = qwe[y==2]
z1 = qwe[y==3]
p0 = p0.apply(lambda x :x.sort_values().values)
p0.drop(p0.index[[0,2]], inplace=True)
p0.drop(p0.index[[-3,-1]],inplace=True)

p1 = p1.apply(lambda x :x.sort_values().values)
p1.drop(p1.index[[0,2]], inplace=True)
p1.drop(p1.index[[-3,-1]],inplace=True)

N1 = N1.apply(lambda x :x.sort_values().values)
N1.drop(N1.index[[0,2]], inplace=True)
N1.drop(N1.index[[-3,-1]],inplace=True)

z1 = z1.apply(lambda x :x.sort_values().values)
z1.drop(z1.index[[0,2]], inplace=True)
z1.drop(z1.index[[-3,-1]],inplace=True)
df = pd.concat([p1,N1,p0,z1],axis = 0)
y = df['y']
X = df.iloc[:,:-1]

#PCA
colors = ['navy', 'turquoise','darkorange','yellow']
target_names = ['P+','N-','P-','Z+']
li = list(range(4210))
wave_lengths = range(4,21)
for wave_length in wave_lengths:    
    start , end = 0,wave_length
    while end < 4210 :
        test = X.iloc[:, start:end]
        left = test.columns[0]
        right = test.columns[-1]
        pca = PCA(n_components= 3)
        x_pca = pca.fit_transform(test)
        x_pca[x_pca >2] = 2
        x_pca[x_pca <-20] = -20
        for pc_number in ([0,1],[1,2],[0,2]):    
            plt.figure(figsize=(20,10))
            for color,i,target_name in zip(colors, [0,1,2,3],target_names):
                plt.scatter(x_pca[y==i,pc_number[0]], x_pca[y==i,pc_number[1]], alpha=.8, color = color,
                           label = target_name)
            plt.legend(loc = 'best', shadow = False, scatterpoints = 1)
            plt.title('%s to %s PCA of raman' %(right,left))
            plt.xlabel('PC%s is %.2f%%' %(pc_number[0],pca.explained_variance_ratio_[pc_number[0]] * 100))
            plt.ylabel('PC%s is %.2f%%'  %(pc_number[1],pca.explained_variance_ratio_[pc_number[1]]  * 100))
            plt.savefig("./picture/%s to %s.png" %(right,left))
            #plt.show()
            start += wave_length
            end += wave_length
