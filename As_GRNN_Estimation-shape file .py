#!/usr/bin/env python
# coding: utf-8

# In[1]:


import pandas as pd 
import numpy as np 
import matplotlib.pyplot as plt
from collections import Counter
import imageio as imageio
from spectral import*
import spectral.io.envi as envi
from collections import Counter
from sklearn.metrics import classification_report


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
# plt.figure(figsize=(10,5))
for i in range(0,33):
    y1=CR[:,i+1]
    first_derivative[:,i] = np.diff(y1,n=1)
fd1=np.zeros((185,34))
fd1[:,0]= wavelen[1:]
fd1[:,1:]=first_derivative

b1=[466,1269,1536,1656,1676,1986,2106]
# b1=[869,1046,1096,1136,1236,1266,1486,1556,1716,1786,2086,2146]
b1=np.array(b1)

FD_resample= np.zeros((len(b1), first_derivative.shape[1]))
a=[]
for i in range(0, len(b1)):
    for j in range (0,185):
        if fd1[j,0] == b1[i] :
            
            FD_resample[i,:]= fd1[j,1:]


# In[4]:


from sklearn.model_selection import  train_test_split
X=FD_resample.T
y=normalized
X_train,X_test , y_train, y_test = train_test_split(X,y , test_size = 0.3, random_state =11)


# In[5]:


from pyGRNN import GRNN
from sklearn.metrics import mean_squared_error as MSE
from sklearn.model_selection import  GridSearchCV


# In[6]:


IGRNN = GRNN()
params_IGRNN = {'kernel':["RBF"],
                'sigma' : list(np.arange(0.1, 4, 0.1)),
                'calibration' : ['warm_start']
                 }
grid_IGRNN = GridSearchCV(estimator=IGRNN,
                          param_grid=params_IGRNN,
                          scoring='neg_mean_squared_error',
                          cv=10,
                          verbose=0
                          )
grid_IGRNN.fit(X_train, y_train.ravel())
best_model = grid_IGRNN.best_estimator_
y_pred = best_model.predict(X_test)
mse_IGRNN = MSE(y_test, y_pred)

mse_IGRNN


# In[7]:


img=open_image('C:/Users/kimmn/Desktop/uni-done/python programming/hyperion/smoothed_cr4.hdr')
img1=np.array(img.open_memmap())
img1[np.isnan(img1)] = 0

b=[467.519989,1265.560059,1537.920044,1659.0,1679.199951,1991.959961,2002.060059]


# In[8]:


img17=np.zeros((702, 428, len(b)))
for i in range (0,161):
    for j in range (0,len(b)):
        if img.bands.centers[i] == b[j]:
#             print(b[j])
            img17[:,:,j]= img1[:,:428,i]
    
from tqdm import tqdm
def extract_pixels(dataset):
    df = pd.DataFrame()
    for i in tqdm(range(dataset.shape[2])):
        df = pd.concat([df, pd.DataFrame(dataset[:, :, i].ravel())], axis=1)
        df = pd.concat([df], axis=1)
    df.columns = [f'band-{i}' for i in range(1, 1+dataset.shape[2])]
    return df


dfimg17 = extract_pixels(img17)

dfimg1=np.array(dfimg17)


# In[9]:


df = dfimg17[dfimg17 != 0]
dfimg2=np.array(df)
dfimg2[np.isnan(dfimg2)] =0


# In[10]:


img_pred = best_model.predict(dfimg2)
np.unique(img_pred ,return_counts=True)


# In[11]:


img_f=np.reshape(img_pred,(702, 428, 1))
imshow(img_f, figsize=(5,10))


# In[12]:


np.unique(img_f , return_counts=True)


# In[13]:


img_pred1 = best_model.predict(dfimg2)
img_f2=np.reshape(img_pred1,(702, 428, 1))


# In[14]:


img_f2[0,0,0]


# In[15]:


pd.DataFrame(img_f2[:,:,0])


# In[16]:


im2=img_f2
for i in range (0,702):
    for j in range(0,428):
        if img_f2[i,j] == 0.09512459703521106  : 
            im2[i,j] = 0.078


# In[17]:


np.unique(im2)


# In[18]:


for i in range (0,702):
    for j in range(0,428):
        im2[i,j] = (im2[i,j]*(max(x)-min(x)))+min(x)


# In[19]:


np.unique(im2)


# In[20]:


plt.figure(figsize=(5,10))
plt.imshow(im2 , cmap='tab20c')
plt.colorbar(fraction=0.07, pad=0.04, label="ppm")
plt.title('Arsenic in soil', fontsize=20)


# In[21]:


# plt.figure(figsize=(10,20))
# plt.imshow(im2 ,cmap='tab20c')
# plt.scatter(shape.x,shape.y)
# plt.imshow(shape)
# plt.colorbar(fraction=0.07, pad=0.04)
# plt.title('Arsenic in soil', fontsize=20)


# In[22]:


import fiona
shapef= fiona.open("C:/Users/kimmn/Desktop/uni-done/arc_epi/patients_hyperion_.shp")


# In[23]:


import shapefile

shape = shapefile.Reader("C:/Users/kimmn/Desktop/uni-done/arc_epi/patients_hyperion_.shp")


# In[24]:


import geopandas as gpd
shape=gpd.read_file("C:/Users/kimmn/Desktop/uni-done/arc_epi/patients_hyperion_.shp")
shape.plot()


# In[25]:


shape


# In[26]:


x=np.array(shape.x)
y=np.array(shape.y)
plt.scatter(x,y)


# In[27]:


plt.figure(figsize=(5,10))
plt.imshow(im2 , cmap='tab20c')
plt.colorbar(fraction=0.07, pad=0.04, label="ppm")
plt.title('Arsenic in soil', fontsize=20)
plt.scatter(x,y )

