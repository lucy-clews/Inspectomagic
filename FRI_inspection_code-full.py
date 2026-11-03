#!/usr/bin/env python
# coding: utf-8

# In[1]:


print('Starting Lucys amazing script :D ')
from astropy.io import fits
import astropy.io.fits as pyfits
from astropy import table
from astropy.table import Table
from matplotlib.colors import LogNorm
import matplotlib


from astropy.visualization import (imshow_norm,AsymmetricPercentileInterval,SqrtStretch,MinMaxInterval)
import matplotlib.image as mpimg
from astropy.utils.data import get_pkg_data_filename
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord  # High-level coordinates
from astropy.coordinates import ICRS, Galactic, FK4, FK5  # Low-level frames
from astropy.coordinates import Angle, Latitude, Longitude  # Angles
import astropy.units as u

import time as time
from numpy import errstate

import matplotlib.pyplot as plt
import matplotlib.image as mpimg

import matplotlib.pyplot as plt

import numpy as np

import sys
import pandas as pd
import glob as glob
import math
import statistics
#import seaborn as sns
import math
import os
import warnings
warnings.filterwarnings('ignore')


print('Python modules loaded')

    
def contour_percentiles(fits_cutout_path, ax, **kwargs):
    
    data = fits.getdata(fits_cutout_path, ext=0) * 1000
    data = data.flatten()
    data = data[~np.isnan(data)]
    percentiles=np.linspace(85,99,10)
    levels = np.percentile(data, percentiles)

    ax.contour(
        data.reshape(fits.getdata(fits_cutout_path,ext=0).shape),
        levels=levels,
        colors='white', linewidths=0.8,alpha=0.5
    )

# In[2]:


def read_csv(filepath):
    
    data = pd.read_csv(filepath, header=0, index_col=0)
    
    return data


def classification(source, acceptable_response, hostxdeg,hostydeg):
    fig = plt.figure(figsize=(6,6))
    fits_cutout_path=glob.glob('/beegfs/lofar/mjh/rgz/*/ridgeline/*/fits_cutouts/%s-cutout.fits' %source)
    hdu = fits.open(fits_cutout_path[0])[0]
    wcs = WCS(hdu.header)
    coords= SkyCoord(hostxdeg, hostydeg, unit='deg')
    hostx,hosty=SkyCoord.to_pixel(coords, wcs)
    ax=fig.add_subplot(1,1,1)
    fits_image=fits.getdata(fits_cutout_path[0],ext=0)*1000
    im,norm=imshow_norm(fits_image, ax, origin='lower',interpolation=None,interval=AsymmetricPercentileInterval(70,100), stretch=SqrtStretch())
    plt.title(source)
    plt.plot(hostx,hosty, 'x', color='r',alpha=0.5)
    cbar=fig.colorbar(im)
    cbar.set_label('mJy/beam', size=12)
    plt.axis('off')    
 
    plt.show(block=False)
    plt.pause(0.5)
 
    description = input("Enter subclass ('a' for lobed/'b' for tailed/'c' for transitional/'d' for half and half/'e' for unclassifiable): ")
    while description not in acceptable_response:
        print("please input one of: 'a', 'b', 'c', 'd','e'")
        description = input("Enter subclass ('a' for lobed/'b' for tailed/'c' for transitional/'d' for half and half/'e' for unclassifiable: ")
  
    plt.close()
    
    if description=='a':
        description='lobed'
    if description=='b':
        description='tailed'
    if description=='c':
        description='transitional'
    if description=='d':
        description='lobe-tail'
    if description=='e':
        description='unclassifiable'
    
    
    return description


acceptable_response =['a', 'b', 'c', 'd','e']

user = input("Please input your initials: ")
save_dir = input("Please enter your working directory path: ")
save_path = f"{save_dir}/{user}_FRI_subclass_env3.0.csv"




if os.path.exists(save_path):
    print("Previous classification file found,  resuming classification...")
    merged_cat = pd.read_csv(save_path, index_col=0)
    remaining = merged_cat[merged_cat['subclass'].isna()]

    print(f"{len(remaining)} sources remaining to classify.")
else:
    print("No previous file found  starting fresh.")
    print('Results will be saved to:', save_path)
    #read in catalogues
    
    with fits.open('/beegfs/lofar/jcroston/surveys/new_envs/lotss_dr2_envs_v0.91.fits') as hdu: # Judith's env catalogue 1st version
        env_catalogue = table.Table(hdu[1].data)
    env_catalogue_old = env_catalogue.to_pandas()
    
    with fits.open('/beegfs/lofar/jcroston/surveys/new_envs/lotss_DR2_AGN_optical_envs_v3.0.fits') as hdu: #Judith's env catalogue latest version
        env_catalogue = table.Table(hdu[1].data)
    env_catalogue_new = env_catalogue.to_pandas()
    print('------------------------------')
    print('Environments catalogue read in, it contains',len(env_catalogue_new),'sources.')
    env_catalogue_new =env_catalogue_new [ (env_catalogue_new['flag_env']==1)&(env_catalogue_new['flag_cen']==1)]
    print('Applying env cat quality flags: \nflag_env=1 (optical richness +ve and L500 frac uncert<1)  \nflag_cen=1 (BCG centre within 250 kpc)')
    print('After quality flags', len(env_catalogue_new), 'sources remain in env catalogue.')
    
    
    df = pd.read_csv('/beegfs/lofar/lclews/DR2_morph_work/paper_1/Data/classifications_v1.csv', low_memory=False)
   
    tmp_cat = df[(df['parent_label']=='FRI')|(df['parent_label']=='bent')]
    FRI_cat = tmp_cat[(tmp_cat['quality']==1)&(tmp_cat['z_best']<=0.8)&(tmp_cat['z_best']>=0.2)]
    print('------------------------------')
    print('FRI catalogue read in')
    
    print('Finding sources which are in both')
    
    cols_old = env_catalogue_old.columns.difference(FRI_cat.columns)
    combined_cat_old = pd.merge(FRI_cat, env_catalogue_old[cols_old], how='inner', left_on='Source_Name_1', right_on='Source_Name')
    cols_new= env_catalogue_new.columns.difference(FRI_cat.columns)
    combined_cat_new = pd.merge(FRI_cat, env_catalogue_new[cols_new], how='inner', left_on='Source_Name_1', right_on='Source_Name')
    
    print('catalogues merged')
    
    #print( 'found',len(combined_cat_old),' sources 0.2<=z<=0.8 in Clews et al 2025 FRIs and old env cat')
    print( 'found',len(combined_cat_new),'Q1 FRIs 0.2<=z<=0.8 in both Clews et al 2025 and new, filtered  env cat.')

    
    #find sources only in NEW catalogue
    env_diff = env_catalogue_new.merge(
        env_catalogue_old[['Source_Name']].drop_duplicates(),
        on='Source_Name',
        how='left',
        indicator=True
    )
    
    env_new_only = env_diff[env_diff['_merge'] == 'left_only'] \
                        .drop(columns='_merge')
    
    
    
    merged_cat = pd.merge(FRI_cat,env_new_only[cols_new],how='inner',left_on='Source_Name_1',right_on='Source_Name')

    if 'subclass' not in merged_cat.columns:
        merged_cat['subclass'] = np.nan
    merged_cat.to_csv(save_path)

    print( "found ", len(merged_cat),"new sources to classify that weren't in old env cat")
    remaining=merged_cat
    print('----------------------------') 
    
print('starting visual inspection of %s sources' %len(remaining))
    
for index, row in remaining.iterrows():
    try:
        source = row['Source_Name']
        print('classifying source:', source)
        hostxdeg=row['optRA']
        hostydeg=row['optDec']
       
        description=classification(source, acceptable_response, hostxdeg,hostydeg)
     
        merged_cat.loc[index, 'subclass'] = description
        merged_cat.to_csv(save_path)
        
    except KeyboardInterrupt:
        print("\n Interrupted. Saving progress...")
        merged_cat.to_csv(save_path)
        done = merged_cat['subclass'].notna().sum()
        total = len(merged_cat)
        print(f"Progress: {done}/{total} sources completed "f"({100*done/total:.1f}%)")
        break
      

print('finished')




