```
import  pandas as pd
import numpy as np
def process(df):
    df['group'] = df['Well'].str.extract(r'([A-Z])') # 将Well 分组
    df['days'] = df['Sample'].str.extract(r'([0-9])')
    df['cq_mean'] = df.groupby(['group','Target'])['Cq'].transform( lambda x : x.mean()) # 按Well 和Target分组求出平均值
    gp = df.groupby('group',as_index=False) # 按A,B,C..分组

    gp_a = gp.apply( lambda x: x.iloc[0:3])
    gp_b = gp.apply( lambda x : x.iloc[3:6])
    gp_c = gp.apply( lambda x : x.iloc[6:9])
    gp_actin = gp.apply( lambda x : x.iloc[9:12])
    for i in [gp_a,gp_b,gp_c]:
        i['actin_mean-cq_mean'] = gp_actin['cq_mean'].values - i['cq_mean'].values
        i['2^(actin_mean-cq_mean)'] = 2**i['actin_mean-cq_mean']
        i['each-actin-cq'] = gp_actin['Cq'].values - i['Cq'].values
        i['each-2^(actin-cq)'] = 2 **i['each-actin-cq']
    result = pd.concat( [gp_a,gp_b,gp_c])
    results = result.sort_values('Well')
    results['std'] = results.groupby(['group', 'Target'])['each-actin-cq'].transform( lambda x : x.std())
    data = results[['Well','Fluor','Target','days','cq_mean','actin_mean-cq_mean','2^(actin_mean-cq_mean)','each-actin-cq','each-2^(actin-cq)','std']]
    '''
    writer = pd.ExcelWriter('.//1.xls')
    data.to_excel(writer)
    writer.save()
    '''
    return data
```
`df= pd.read_excel('E.gracilis(2)_2018-12-30 14-49-50_CT031429 -  Quantification Cq Results.xls')`     
`data = process(df);
writer = pd.ExcelWriter('.//2018-12-30-3.xls')
data.to_excel(writer)
writer.save()`
