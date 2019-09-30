#!/usr/bin/env python
# coding: utf-8

# # To generate the data using in Figure S3 of the supplement

# In[12]:


import scipy.io
from scipy import stats, optimize, interpolate
import pandas as pd
import numpy as np
from collections import defaultdict
import math
import matplotlib.pyplot as plt
import matplotlib.cm as cm
get_ipython().run_line_magic('matplotlib', 'inline')


# ## import the data from mat format and generate a dataframe with it


# In[29]:


names=['d1','d2','d3','d4','d5','d6','d7']
meaning=['0.05% Gal','0.2% Gal','2% Gal','2% Glucose','4% Gal','Cdc42 promoter','NO-GFP']
reads=defaultdict(dict)
i=0
for each in names:
    r= scipy.io.loadmat(each + '.mat')
    reads[each]['data'] =r[each]
    reads[each]['log10']=np.log10(r[each])
    reads[each]['geomean_log10']=scipy.stats.mstats.gmean(np.log10(r[each]), axis=0)
    reads[each]['geomean']=scipy.stats.mstats.gmean((r[each]), axis=0)
    reads[each]['std_log10']=np.std(np.log10(r[each]))
    reads[each]['Names']=meaning[i]
    reads[each]['Length of data']=len(r[each])
    i=i+1
pd_reads=pd.DataFrame(reads).T


# In[30]:


pd_reads


# In[31]:


pd_reads.to_excel('data-histograms.xlsx')


# ## Figure S3 histograms

# In[33]:


fig, axes = plt.subplots(1,1, figsize=(10,5), dpi=100, sharex=True, sharey=True)
kwargs = dict(alpha=0.5, bins=20, density=True, stacked=True)
axis_font = {'fontname':'Arial', 'size':'12'}
colors = cm.rainbow(np.linspace(0, 1, len(pd_reads['Names'])))
keys=pd_reads['Names']   
i=0
for y, c in zip(keys, colors):
    plt.hist(pd_reads['log10'][i],bins=20,density=True,label=pd_reads['Names'][i],stacked=True,alpha=0.5,color=c)
    i=i+1


#plt.hist(pd_reads['log10'],**kwargs,label=pd_reads['Names'])
plt.legend(prop={'size': 10},loc='upper left')
plt.xlabel('log(10)_Intensities',**axis_font)
plt.ylabel('Normalized frequency',**axis_font)


# In[ ]:





# ## Bar plot of the geometric mean and standard deviations (inset of the Fig S3)

# In[5]:


fig, axes = plt.subplots(1,1, figsize=(10,5), dpi=100, sharex=True, sharey=True)
kwargs = dict(alpha=0.5, bins=20, density=True, stacked=True)
axis_font = {'fontname':'Arial', 'size':'12'}
plt.bar(x=pd_reads['Names'],height=pd_reads['geomean_log10'],yerr=pd_reads['std_log10'],align='center', alpha=0.5, ecolor='black', capsize=10)
plt.xlabel('log(10)_Intensities',**axis_font)
plt.ylabel('Geometric mean of Log10(Intensities)',**axis_font)


# In[ ]:





# In[ ]:




