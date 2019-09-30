print("hello world")

import scipy.io
import scipy
from scipy import stats, optimize, interpolate
import pandas as pd
import numpy as np
from collections import defaultdict
import math
import matplotlib.pyplot as plt

names=['data_point_1','data_point_2','data_point_3','data_point_4','data_point_5','data_point_6']
variable_names=['d1','d2','d3','d4','d5','d6']

time_points=[1,3,6,16,18,21]

reads=defaultdict(dict)
i=0
for each in names:
    r= scipy.io.loadmat("C:\\Users\\linigodelacruz\\Documents\\PhD_2018\\Work with Fridtjof paper\\Calculations-supplement\\Figure-S2-mean-intensity- Gal-Cdc42-GFP-promoter\\" + each + '.mat')
    reads[each]['data'] =r[variable_names[i]]
    reads[each]['log10']=np.log10(r[variable_names[i]])
    reads[each]['time-points']=time_points[i]
    reads[each]['geomean_log10']=scipy.stats.mstats.gmean(np.log10(r[variable_names[i]]), axis=0)
    reads[each]['std_log10']=np.std(np.log10(r[variable_names[i]]))
    reads[each]['Length of data']=len(r[variable_names[i]])
    reads[each]['Names']=(variable_names[i])
    i=i+1
pd_reads=pd.DataFrame(reads).T



# ## Lne plot with error bar of figure S2 



fig, axes = plt.subplots(2,1, figsize=(10,5), dpi=100, sharex=True, sharey=True)
kwargs = dict(alpha=0.5, bins=20, density=True, stacked=True)
axis_font = {'fontname':'Arial', 'size':'10'}
plt.subplot(2,1,1)
plt.errorbar(x=pd_reads['time-points'],y=np.power(pd_reads['geomean_log10'],10),yerr=np.power(pd_reads['std_log10'],10), alpha=0.5, ecolor='black', capsize=10)
plt.grid(True)

plt.ylabel('Geometric mean of (Intensities)',**axis_font)
plt.subplot(2,1,2)
plt.errorbar(x=pd_reads['time-points'],y=pd_reads['geomean_log10'],yerr=pd_reads['std_log10'], alpha=0.5, ecolor='black', capsize=10)
plt.grid(True)
plt.ylabel('Geometric mean of Log(Intensities)',**axis_font)
plt.xlabel('Time (hours)')
plt.show()
