

# In[1]:


import scipy.io
import scipy
from scipy import stats, optimize, interpolate
import pandas as pd
import numpy as np
from collections import defaultdict
import math
import matplotlib.pyplot as plt
import matplotlib.cm as cm

# ## Creating the data arrays from the estimated and measured quantities

# In[5]:


names=['d1','d2','d3','d5','d6','d7','d4']
meaning=['0.05% Gal','0.2% Gal','2% Gal','4% Gal','Cdc42 promoter','NO-GFP','2% Glucose']

N0=8690
reads=defaultdict(dict)
i=0
# measured quantities
for each in names:
    r= scipy.io.loadmat("C:\\Users\\linigodelacruz\\Documents\\PhD_2018\\Work with Fridtjof paper\\Calculations-supplement\\Figure-S6\\" + each + '.mat')
    reads[each]['data'] =r[each]
    reads[each]['log10_data'] =np.log10(r[each])
    reads[each]['data_std'] =np.std(r[each])
    reads[each]['geomean']=scipy.stats.mstats.gmean((r[each]), axis=0)
    reads[each]['Names']=meaning[i]
    #reads[each]['Length of data']=len(r[each])
    reads[each]['geomean_log10']=scipy.stats.mstats.gmean(np.log10(r[each]), axis=0)
    
    i=i+1
# estimated quantities
for each in names:
    reads[each]['cdc42_copy'] = N0*(reads[each]['geomean']-reads['d4']['geomean'])/(reads['d6']['geomean']-reads['d4']['geomean'])
    reads[each]['cdc42_copy/N0'] = (reads[each]['geomean']-reads['d4']['geomean'])/(reads['d6']['geomean']-reads['d4']['geomean'])
    reads[each]['cdc42_hist/N0'] = (reads[each]['data']-reads['d4']['geomean'])/(reads['d6']['geomean']-reads['d4']['geomean'])
    reads[each]['cdc42_hist-std'] = np.std(reads[each]['data']-reads['d4']['geomean'])/(reads['d6']['geomean']-reads['d4']['geomean'])
    reads[each]['cdc42_hist-ste'] = scipy.stats.sem(reads[each]['data']-reads['d4']['geomean'])/(reads['d6']['geomean']-reads['d4']['geomean'])
    reads[each]['cdc42_copy_8h'] = N0*(0.3*reads[each]['geomean']-reads['d4']['geomean'])/(reads['d6']['geomean']-reads['d4']['geomean'])



# ## Histograms of the cdc42 molecules

# In[7]:


pd_reads=pd.DataFrame(reads).T
fig, axes = plt.subplots(1,1, figsize=(10,5), dpi=100, sharex=True, sharey=True)
kwargs = dict(alpha=0.5, bins=20, density=True, stacked=True)
axis_font = {'fontname':'Arial', 'size':'12'}
colors = cm.Blues(np.linspace(0.5, 1.2, len(pd_reads['Names'])))
keys=pd_reads['Names'][1:4] 
i=1

for y, c in zip(keys, colors):
    plt.hist(pd_reads['cdc42_hist/N0'][i],bins=50,density=False,label=pd_reads['Names'][i],stacked=True,alpha=0.7,color=c)
    i=i+1


#plt.hist(pd_reads['log10'],**kwargs,label=pd_reads['Names'])
#plt.xlim(-0.5,17)
plt.legend(prop={'size': 10},loc='upper right')
plt.xlabel('Normalized cdc42 number by the reported normal number (8690) in WT',**axis_font)
plt.ylabel('Normalized frequency',**axis_font)
#plt.xscale('log')

#plt.savefig("cdc42-by-intensities.png",dpi=300,format='png',transparent=True)
#plt.show()

## Boxplots of cdc42 distributions 

cdc42_dist=[]
for i in np.arange(0,4):
    cdc42_dist.append(pd_reads['cdc42_hist/N0'][i])
cdc42_copy=[]
for i in np.arange(0,4):
    cdc42_copy.append(pd_reads['cdc42_copy/N0'][i])
cdc42_names=[]
for i in np.arange(0,4):
    cdc42_names.append(pd_reads['Names'][i])


fig, axes = plt.subplots(1,1, figsize=(10,5), dpi=100, sharex=True, sharey=True)
kwargs = dict(alpha=0.5, bins=20, density=True, stacked=True)
axis_font = {'fontname':'Arial', 'size':'12'}
colors = cm.Blues(np.linspace(0.5, 1.2, len(pd_reads['Names'])))
keys=pd_reads['Names'] 
i=0

b=plt.boxplot(cdc42_dist)
p=plt.plot([1,2,3,4],cdc42_copy,color='w', marker='*', markeredgecolor='b',MarkerSize=8)

#plt.hist(pd_reads['log10'],**kwargs,label=pd_reads['Names'])
#plt.xlim(-0.5,17)
plt.ylabel('cdc42/N0',**axis_font)

fig.text(0.715,0.8, ' * Geometric Mean', color='blue', weight='roman',
         size='x-small')

#plt.xscale('log')

plt.xticks([1,2,3,4], cdc42_names)
plt.show()



## Histograms of log of the Cdc42
fig, axes = plt.subplots(1,1, figsize=(10,5), dpi=100, sharex=True, sharey=True)
kwargs = dict(alpha=0.5, bins=20, density=True, stacked=True)
axis_font = {'fontname':'Arial', 'size':'12'}
colors = cm.Blues(np.linspace(0.2, 1.2, len(pd_reads['Names'])))
keys=pd_reads['Names']   
i=0

for y, c in zip(keys, colors):
    plt.hist(pd_reads['log10_data'][i],bins=50,density=True,label=pd_reads['Names'][i],stacked=True,alpha=0.7,color=c)
    i=i+1


#plt.hist(pd_reads['log10'],**kwargs,label=pd_reads['Names'])
#plt.ylim(0,4000)
plt.legend(prop={'size': 10},loc='upper left')
plt.xlabel('Raw data intensities',**axis_font)
plt.ylabel('Normalized frequency',**axis_font)
plt.xlim(0.5,3.7)
plt.savefig("Figure-S3.png",dpi=300,format='png',transparent=True)


# ## cdc42 copy number calculation

# In[78]:


names=['data_doublings','std_doublings'] # estimation from the population growth data
# data_doublings=[doubling @0%,doubling @0.06%,doubling @0.1%]

data_doubling_2_percent=[3/0.006,1/0.006,2/0.009,1/0.009,1/0.011,1/0.011]
meaning=['data_8h_bem1d','data_8h-24h_bem1d','data_8h_bem1dbem3d','data_8h-24h_bem1dbem3d','data_8h_wt','data_8h-24h_wt']

strains_8h=['bem1d-8h','bem1dbem3d-8h','WT+cdc42-8h']
strains_24h=['bem1d-24h','bem1dbem3d-24h','WT+cdc42-24h']

calculation=defaultdict(dict)
pos=[0,2,4] # the positions of 8h data points for the doubling times data

    
r= scipy.io.loadmat("C:\\Users\\linigodelacruz\\Documents\\PhD_2018\\Work with Fridtjof paper\\Calculations-supplement\\\Figure-S6\\" + names[0] + '.mat')#dooubling times in min
r_std= scipy.io.loadmat("C:\\Users\\linigodelacruz\\Documents\\PhD_2018\\Work with Fridtjof paper\\Calculations-supplement\\\Figure-S6\\"+ names[1] + '.mat')#std doubling times in mins

# I have to fix the doubling times as a combination of the growth rate in 0% ,0.06% and 0.1%and in 2% , so far it is just for the
#specific concentrations
number_of_doublings=np.zeros((3,6))
for i in pos:
    for j in np.arange(0,3):
        if r[names[0]][:, i][j]==0:
            number_of_doublings[j,i]=8*60/data_doubling_2_percent[i]*0.5
        else:
             number_of_doublings[j,i]=np.mean([8*60/r[names[0]][:, i][j],8*60/data_doubling_2_percent[i]])



for each,i in zip(strains_8h,pos): 

    calculation[each]['Names'] = each
    calculation[each]['# of doublings'] = number_of_doublings[:,i]
    calculation[each]['std # of doublings'] = r_std[names[1]][:, i]/60
    #calculation[each]['0%Gal-std'] = reads['d3']['cdc42_hist-std']/(np.power(2,r_std[names[1]][:, i][0]/60)) #converting to hours
    calculation[each]['0%Gal'] = 1/2*reads['d3']['cdc42_copy_8h']/(np.power(2,number_of_doublings[:, i][0])) #converting to hours
    calculation[each]['0.06%Gal'] =1/2*reads['d3']['cdc42_copy_8h']/(np.power(2,number_of_doublings[:, i][1]))+reads['d1']['cdc42_copy_8h']
    calculation[each]['0.1%Gal'] = 1/2*reads['d3']['cdc42_copy_8h']/(np.power(2,number_of_doublings[:, i][2]))+0.5*reads['d2']['cdc42_copy_8h']

pos_24h=[1,3,5]# the positions of 24h data points for the doubling times data

for i in pos_24h:
    for j in np.arange(0,3):
        if r[names[0]][:, i][j]==0:
            number_of_doublings[j,i]=24*60/data_doubling_2_percent[i]*0.5
        else:
             number_of_doublings[j,i]=np.mean([24*60/r[names[0]][:, i][j],24*60/data_doubling_2_percent[i]])
        
        
for each,i in zip(strains_24h,pos_24h): 
    
    calculation[each]['Names'] = each
    calculation[each]['# of doublings'] = number_of_doublings[:,i]
    calculation[each]['std # of doublings'] = r_std[names[1]][:, i]/60
    calculation[each]['0%Gal'] = 1/8*reads['d3']['cdc42_copy']/(np.power(2,number_of_doublings[:, i][0]))
    calculation[each]['0.06%Gal'] = 1/8*reads['d3']['cdc42_copy']/(np.power(2,number_of_doublings[:, i][1]))+reads['d1']['cdc42_copy']
    calculation[each]['0.1%Gal'] = 1/8*reads['d3']['cdc42_copy']/(np.power(2,number_of_doublings[:, i][2]))+0.5*reads['d2']['cdc42_copy']

strains_wt=['WT-8h','WT-24h']  
for each in strains_wt:
    calculation[each]['0%Gal']=8690
    calculation[each]['0.06%Gal']=8690
    calculation[each]['0.1%Gal']=8690
    calculation[each]['Names'] = each


# In[34]:


[3/0.006,1/0.006,2/0.009,1/0.009,1/0.011,1/0.011]


# In[79]:




pd_calculation=pd.DataFrame(calculation).T
pd_calculation.T


# In[54]:


#measured quantities 
names=['bem1d + Galpr-CDC42 @8h','bem1bem3 deleted + Galpr-CDC42 @8h','WT + Galpr-CDC42 @8h','WT @8h','bem1d + Galpr-CDC42 @24h',
      'bem1bem3 deleted + Galpr-CDC42 @24h','WT + Galpr-CDC42 @24h','WT @24h']
cell_radii=np.array([[5.825, 4.240, 4.125],[5.25, 3.995 , 3.915],[3.185,3.25,3.2],[3.29,3.13, 3.145],
                    [7.36,4.45,4.22],[7.5,3.6,3.97],[6.12,3.04,3.35],[2.85,2.94,2.96]]) # measured quantity



cell_radii_std=np.array([[0.601040764,0.22627417,0.318198052],[0.494974747,0.417193001,0.16263456],
                         [0.06363961,0.070710678,0.141421356],[0.296984848,0,0.06363961],[0.93,1.33,1.68],
                         [0.75,1.11,1.85],[0.93,1.29,0.75],[0.43,0.42,0.83]])


# In[55]:


quantities_0=defaultdict(dict)
i=0
names=['b1d-8h','b1-3d-8h','WT+42-8h','WT-8h','b1d-24h',
      'b1-3d-24h','WT+42-24h','WT-24h']
strains=['bem1d-8h','bem1dbem3d-8h','WT+cdc42-8h','WT-8h','bem1d-24h','bem1dbem3d-24h','WT+cdc42-24h','WT-24h']


# a forloop per galactose concentration 
for strains_name in strains:
        quantities_0[strains_name]['Name']=names[i]
        quantities_0[strains_name]['cdc42_copy_number']=calculation[strains_name]['0%Gal']                                                                     
        quantities_0[strains_name]['cell_radius(um)']=cell_radii[i,:][0]
        quantities_0[strains_name]['cell_radius-std']=cell_radii_std[i,:][0]
        quantities_0[strains_name]['cell volume']=4/3*math.pi*np.power(cell_radii[i,:][0],3)
        quantities_0[strains_name]['cdc42 density']=calculation[strains_name]['0%Gal']/quantities_0[strains_name]['cell volume']
        i=i+1
quantities_1=defaultdict(dict) #0.06%
i=0
for strains_name in strains:
        quantities_1[strains_name]['Name']=names[i]
        quantities_1[strains_name]['cdc42_copy_number']=calculation[strains_name]['0.06%Gal']                                                                     
        quantities_1[strains_name]['cell_radius(um)']=cell_radii[i,:][1]
        quantities_1[strains_name]['cell_radius-std']=cell_radii_std[i,:][1]
        quantities_1[strains_name]['cell volume']=4/3*math.pi*np.power(cell_radii[i,:][1],3)
        quantities_1[strains_name]['cdc42 density']=calculation[strains_name]['0.06%Gal']/quantities_1[strains_name]['cell volume']
        i=i+1
        
quantities_2=defaultdict(dict) #0.06%
i=0
for strains_name in strains:
        quantities_2[strains_name]['Name']=names[i]
        quantities_2[strains_name]['cdc42_copy_number']=calculation[strains_name]['0.1%Gal']                                                                     
        quantities_2[strains_name]['cell_radius(um)']=cell_radii[i,:][2]
        quantities_2[strains_name]['cell_radius-std']=cell_radii_std[i,:][2]
        quantities_2[strains_name]['cell volume']=4/3*math.pi*np.power(cell_radii[i,:][2],3)
        quantities_2[strains_name]['cdc42 density']=calculation[strains_name]['0.1%Gal']/quantities_1[strains_name]['cell volume']
        i=i+1


# In[56]:


keys=['0%-Gal','0.06%-Gal','0.1%-Gal']
pd_quantities_0=pd.DataFrame(quantities_0).T
pd_quantities_1=pd.DataFrame(quantities_1).T
pd_quantities_2=pd.DataFrame(quantities_2).T

frames = [pd_quantities_0,pd_quantities_1,pd_quantities_2]
df_quantities = pd.concat(frames, keys=keys,sort=False).T


df_quantities.T


# In[57]:


pd_reads=df_quantities.T
pd_reads.T[keys].loc['Name']


# In[70]:


pd_reads=df_quantities.T
gal=[['0%-Gal'],['0.06%-Gal'],['0.1%-Gal']]
fig, axes = plt.subplots(1,3, figsize=(10,10), dpi=100, sharex=True, sharey=True)
plt.subplots_adjust(bottom=0, right=2.5, top=0.7)
kwargs = dict(alpha=0.5, bins=20, density=True, stacked=True)
axis_font = {'fontname':'Arial', 'size':'20'}
#colors = cm.rainbow(np.linspace(0, 1, len(pd_reads[0:4])))
colors=['blue','black','lightgreen']

variables=['cdc42 density','cdc42_copy_number','cell volume']
width = 0.5 
j=1
for keys, c in zip(gal,colors): 
    plt.subplot(1,3,j)
    plt.bar(x=pd_reads.T[keys].loc['Name'],height=pd_reads.T[keys].loc[variables[0]],alpha=0.8,width=width,color=c)
    plt.title(keys,**axis_font)
    plt.tick_params(labelsize=14)
    plt.ylabel(variables[0],**axis_font)
    plt.ylim(0,100)
    j=j+1

fig, axes = plt.subplots(1,3, figsize=(10,10), dpi=100, sharex=True, sharey=True)
plt.subplots_adjust(bottom=0, right=2.5, top=0.7)
j=1
for keys, c in zip(gal,colors): 
    plt.subplot(1,3,j)
    plt.bar(x=pd_reads.T[keys].loc['Name'],height=pd_reads.T[keys].loc[variables[1]]/8690,alpha=0.8,width=width,color=c)
    plt.title(keys,**axis_font)
    plt.tick_params(labelsize=14)
    plt.ylabel(variables[1] +'-normalized to the WT',**axis_font)
    plt.ylim(0,1.3)
    j=j+1  
    
fig, axes = plt.subplots(1,3, figsize=(10,10), dpi=100, sharex=True, sharey=True)
plt.subplots_adjust(bottom=0, right=2.5, top=0.7)
j=1
for keys, c in zip(gal,colors): 
    plt.subplot(1,3,j)
    plt.bar(x=pd_reads.T[keys].loc['Name'],height=pd_reads.T[keys].loc[variables[2]],alpha=0.8,width=width,color=c)
    plt.title(keys,**axis_font)
    plt.tick_params(labelsize=14)
    plt.ylabel(variables[2],**axis_font)
    plt.ylim(0,1800)
    j=j+1  


# In[64]:


pd_reads.T[keys].loc[variables[0]].max()


# In[59]:


pd_reads=df_quantities.T
gal=['0%-Gal','0.06%-Gal','0.1%-Gal']
fig, axes = plt.subplots(3,1, figsize=(10,10), dpi=100, sharex=True, sharey=True)
plt.subplots_adjust(bottom=0, right=1, top=1.3)
kwargs = dict(alpha=0.5, bins=20, density=True, stacked=True)
axis_font = {'fontname':'Arial', 'size':'12'}
#colors = cm.rainbow(np.linspace(0, 1, len(pd_reads[0:4])))
colors=['blue','black','lightgreen']

i=0
width = 0.35 
for keys, c in zip(gal,colors): 
    plt.subplot(3,1,1)
    plt.bar(x=pd_reads.T[keys].loc['Name'],height=pd_reads.T[keys].loc['cdc42 density'],alpha=0.3,width=width,color=colors[i],label=keys)
    plt.legend(prop={'size': 10},loc='upper right')
    plt.tick_params(labelsize=8)
    plt.ylabel('Cdc42 density',**axis_font)
    plt.subplot(3,1,2)
    plt.bar(x=pd_reads.T[keys].loc['Name'],height=pd_reads.T[keys].loc['cell volume'],alpha=0.3,width=width,color=colors[i],label=keys)
    plt.ylabel('Cell Volume(um^3)',**axis_font)
    plt.legend(prop={'size': 10},loc='upper right')
    plt.tick_params(labelsize=8)
    plt.subplot(3,1,3)
    plt.bar(x=pd_reads.T[keys].loc['Name'],height=pd_reads.T[keys].loc['cdc42_copy_number'],alpha=0.3,width=width,color=colors[i],label=keys)
    plt.ylabel('Cdc42 copy number',**axis_font)
    plt.legend(prop={'size': 10},loc='upper right')
    plt.tick_params(labelsize=8)
    i=i+1











# In[ ]:





# ## Previous calculation .... not precise

# In[91]:




strains=['bem1d-8h','bem1dbem3d-8h','WT+cdc42-8h','WT-8h', 'bem1d-24h','bem1dbem3d-24h','WT+cdc42-24h','WT-24h']
names=['bem1d + Galpr-CDC42 @8h','bem1bem3 deleted + Galpr-CDC42 @8h','WT + Galpr-CDC42 @8h','WT @8h','bem1d + Galpr-CDC42 @24h',
      'bem1bem3 deleted + Galpr-CDC42 @24h','WT + Galpr-CDC42 @24h','WT @24h']
cdc42_number=np.array([[17403.24309,19488.8430,25512.73588 ],[ 10565.46437,9056.206906,18937.32872],[8232.235602 , 4449.688873,15036.78117],
                       [8690 , 8690,8690],[4768.488607,6958.508065,13958.24344],[2894.937236,2155.662341,13050.67689],[2255.632555,2090.102004,13036.95961],
                       [8690 , 8690, 8690]]) # estimated quantity
#cdc42_number=df.iloc[0,:]       
cell_radii=np.array([[5.825, 4.240, 4.125],[5.25, 3.995 , 3.915],[3.185,3.25,3.2],[3.29,3.13, 3.145],
                    [7.36,4.45,4.22],[7.5,3.6,3.97],[6.12,3.04,3.35],[2.85,2.94,2.96]]) # measured quantity
cell_radii_std=np.array([[0.601040764,0.22627417,0.318198052],[0.494974747,0.417193001,0.16263456],
                         [0.06363961,0.070710678,0.141421356],[0.296984848,0,0.06363961],[0.93,1.33,1.68],
                         [0.75,1.11,1.85],[0.93,1.29,0.75],[0.43,0.42,0.83]])


# In[ ]:





# ## Creating the dataframe of data from the estimated and measured quantities, above defined

# In[92]:


reads=defaultdict(dict)
i=0
indexes=['Gal percentages','cdc42_copy_number','cell_radius(um)','cell_radius-std','cell volume(um^3)','cdc42 density(um^(-3))']
for strains_name in strains:
        reads[strains_name]['Gal percentages']=['0%','0.06%','0.1%']
        reads[strains_name]['Name']=names[i]
        reads[strains_name]['cdc42_copy_number']=cdc42_number[i,:]
        reads[strains_name]['cell_radius(um)']=cell_radii[i,:]
        reads[strains_name]['cell_radius-std']=cell_radii_std[i,:]
        reads[strains_name]['cell volume(um^3)']=4/3*math.pi*np.power(cell_radii[i,:],3)
        reads[strains_name]['cdc42 density(um^(-3))']=cdc42_number[i,:]/reads[strains_name]['cell volume(um^3)']
        i=i+1


# In[182]:


df=pd.DataFrame(reads)
df.T


# ## Save the data as a dataframe to a txt file 

# In[183]:


df.to_excel('data-fig-s6.xlsx')


# In[187]:


df_from_excel=pd.read_excel('data-fig-s6.xlsx',dtype='object')


# In[188]:


df_from_excel.T
#issue: I have a problem when reading the exported dataframe , the format it is not the same as before :(


# ## Selecting the relevant data for the plots

# In[179]:


# cell radius plots


pd_reads=df.T

fig, axes = plt.subplots(2,4, figsize=(10,5), dpi=100, sharex=True, sharey=True)
plt.subplots_adjust(bottom=0.2, right=2, top=1.4)
axis_font = {'fontname':'Arial', 'size':'14','fontweight':'bold'}
gal_percentages=['0%','0.06%','0.1%']
for i in range(8):
    plt.subplot(2,4,i+1)
    plt.bar(x=gal_percentages,height=pd_reads['cell_radius(um)'][i],yerr=pd_reads['cell_radius-std'][i],align='center', alpha=0.5, ecolor='black', color='red',capsize=10)
    plt.ylim(0,9)
    plt.title(pd_reads['Name'][i],**axis_font)
    plt.ylabel('Mean of cell radii',**axis_font)
    


# In[153]:


# cdc42 copy number

pd_reads=df.T

fig, axes = plt.subplots(2,4, figsize=(10,5), dpi=100, sharex=True, sharey=True)
plt.subplots_adjust(bottom=0.2, right=2, top=1.4)
kwargs = dict(alpha=0.5, bins=20, density=True, stacked=True)
axis_font = {'fontname':'Arial', 'size':'14','fontweight':'bold'}

for i in range(8):
    plt.subplot(2,4,i+1)
    plt.bar(x=pd_reads['Gal-percentages'][i],height=pd_reads['cdc42_copy_number'][i], alpha=0.5)
    plt.ylim(0,28000)
    plt.title(pd_reads['name'][i],**axis_font)
    plt.ylabel('cdc42_copy_number',**axis_font)
    
    


# In[92]:


# cell density

pd_reads=df.T
fig, axes = plt.subplots(2,4, figsize=(10,5), dpi=100, sharex=True, sharey=True)
plt.subplots_adjust(bottom=0.2, right=2, top=1.4)
kwargs = dict(alpha=0.5, bins=20, density=True, stacked=True)
axis_font = {'fontname':'Arial', 'size':'14','fontweight':'bold'}

for i in range(8):
    plt.subplot(2,4,i+1)
    plt.bar(x=pd_reads['Gal-percentages'][i],height=pd_reads['cdc42 density(um^(-3))'][i], alpha=0.5,color='green')
    plt.ylim(0,120)
    plt.title(pd_reads['name'][i],**axis_font)
    plt.ylabel('cdc42_density',**axis_font)


# In[ ]:




