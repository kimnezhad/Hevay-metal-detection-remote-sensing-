#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np 
import pandas as pd 
import matplotlib.pyplot as plt
import scipy.stats #for derivatives
from scipy.stats import spearmanr
from scipy import stats
import warnings
# warnings.filterwarnings('ignore')


# In[2]:


CR=np.array(pd.read_excel("CRwithWL.xlsx")) 
name_id=np.array(pd.read_excel("name_and_id_cr.xlsx"))
As_T= np.array(pd.read_excel("As.xlsx"))

wavelen=CR[:,0]
name=name_id[:,1]
As= As_T[:,1]

## normalize AS
normalized=np.zeros((33,1))
x=As
for i in range(0,33):
    normalized[i] = (x[i]-min(x))/(max(x)-min(x))


# In[3]:


first_derivative=np.zeros((185,33))
plt.figure(figsize=(20,10))
for i in range(0,33):
    y1=CR[:,i+1]
    first_derivative[:,i] = np.diff(y1,n=1)
    plt.plot(wavelen[1:],first_derivative[:,i])
    plt.legend(name)
    plt.title("First Derivative CR spectra", fontsize=25)
    plt.xlabel("wavelength(nanometer)",fontsize=20)
    plt.ylabel("reflectance",fontsize=20)
    
df1=first_derivative.T
spectra_df1=pd.DataFrame(np.hstack((normalized,df1)))
corrMatrix_df1 =np.array(spectra_df1.corr(method='spearman')) 
plt.figure(figsize=(20,10))
plt.plot(wavelen[1:],corrMatrix_df1[1:,0])
plt.axhline(y=0, color='r', linestyle='-')
plt.axhline(y=0.35, color='g', linestyle='-')
# plt.axhline(y=-0.2, color='g', linestyle='-')
plt.xlabel("wavlength",fontsize=20)
plt.ylabel("correlation",fontsize=20)
plt.title(" Corelation of First derivative CR spectra with AS content",fontsize=20)
plt.grid()   


# In[4]:


# wavelen[2:].shape


# In[5]:


corr_df_bands=np.zeros((185,2))
corr_df_bands[:,0]=wavelen[1:]
corr_df_bands[:,1]=corrMatrix_df1[1:,0]
pd.DataFrame(corr_df_bands)


# In[6]:


plt.figure(figsize=(20,10))
plt.plot(wavelen[1:],np.abs(corrMatrix_df1[1:,0]))
plt.axhline(y=0.2, color='g', linestyle='-')


# In[7]:


wavelen_cor_As=[]
for i in range(1,185):
    if corr_df_bands[i,1] > 0.35 :
        wavelen_cor_As=np.append(wavelen_cor_As,corr_df_bands[i,0])
        print("the correlation for WL ",corr_df_bands[i,0]," is ",corr_df_bands[i,1]  )
        


# In[8]:


wavelen_cor_As.shape


# In[9]:


ttest_df=stats.ttest_1samp(corrMatrix_df1,0.1)
pval_df=ttest_df.pvalue
pd.DataFrame(ttest_df)
# a=0
# alpha = 0.1

# for i in range(0,185):
#     a=a+1
#     if pval_df[i] < alpha :
#           print(a , 'correlation between these two variables is statistically significant  p=%.3f' % pval_df[i] )
#     else:
#          print(a, 'Samples correlation between these two variables is NOT statistically significant  p=%.3f' % pval_df[i] )


# In[10]:


pval_df2=np.zeros((186,2))
pval_df2[:,0]=wavelen
pval_df2[:,1]= pval_df.T
# pd.DataFrame(pval_df2)


# In[11]:


wavelen_ttest_As=[]
alpha=0.1
for i in range(0,185):
    
    if pval_df2[i,1] < alpha :
        wavelen_ttest_As=np.append(wavelen_ttest_As,pval_df2[i,0])
        
        print('correlation between ' , pval_df2[i,0], 'and As is significant  p=%.3f' % pval_df2[i,1] )
#     else:
#          print('Samples correlation between these two variables is NOT statistically significant  p=%.3f' % pval_df2[i,1] )

