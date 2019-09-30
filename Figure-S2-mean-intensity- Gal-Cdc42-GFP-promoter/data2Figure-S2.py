
# In[1]:


import scipy.io
from scipy import stats, optimize, interpolate
import scipy
import pandas as pd
import numpy as np
from collections import defaultdict
import math
import matplotlib.pyplot as plt



# ## import the data from mat format and generate a dataframe with it

# In[2]:


names=['data_point_1','data_point_2','data_point_3','data_point_4','data_point_5','data_point_6']
variable_names=['d1','d2','d3','d4','d5','d6']

time_points=[1,3,6,16,18,21]

reads=defaultdict(dict)
i=0
for each in names:
    r= scipy.io.loadmat("C:\\Users\\linigodelacruz\\Documents\\PhD_2018\\Work with Fridtjof paper\\Calculations-supplement\\Figure-S2-mean-intensity- Gal-Cdc42-GFP-promoter\\" + each + '.mat')
    reads[each]['data'] =r[variable_names[i]]
    reads[each]['std'] =np.std(r[variable_names[i]])
    reads[each]['ste']=scipy.stats.sem(r[variable_names[i]])
    reads[each]['log10']=np.log10(r[variable_names[i]])
    reads[each]['time-points']=time_points[i]
    reads[each]['geomean_log10']=scipy.stats.mstats.gmean(np.log10(r[variable_names[i]]), axis=0)
    reads[each]['geomean']=scipy.stats.mstats.gmean((r[variable_names[i]]), axis=0)
    reads[each]['mean']=np.mean((r[variable_names[i]]), axis=0)
    reads[each]['log10_geomean']=np.log10(scipy.stats.mstats.gmean((r[variable_names[i]]), axis=0))
    reads[each]['std_log10']=np.std(np.log10(r[variable_names[i]]))
    reads[each]['Length of data']=len(r[variable_names[i]])
    reads[each]['Names']=(variable_names[i])
    i=i+1
pd_reads=pd.DataFrame(reads).T



# In[10]:


fig, axes = plt.subplots(1,1, figsize=(10,5), dpi=100, sharex=True, sharey=True)
kwargs = dict(alpha=0.5, bins=20, density=True, stacked=True)
axis_font = {'fontname':'Arial', 'size':'10'}

plt.errorbar(x=pd_reads['time-points'],y=(pd_reads['geomean']),yerr=pd_reads['std'], alpha=0.5, ecolor='black', capsize=5,fmt='bo')
plt.ylabel('Geometric mean of intensity units')
plt.xlabel('Time (hours)')
plt.grid(True,axis='both')
plt.tight_layout()
plt.show()
#plt.savefig("Figure-S2.png",dpi=300,format='png',transparent=True)


