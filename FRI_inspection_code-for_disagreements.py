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

import matplotlib.pyplot as plta
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
    try:
        ax.contour(
        data.reshape(fits.getdata(fits_cutout_path, ext=0).shape),
        levels=levels,
        colors='white', linewidths=0.8, alpha=0.3
    )
    except:
        print('COUNTOURS FAILED')
  

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
    contour_percentiles( fits_cutout_path[0], ax)
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
old_col= f'subclass_{user}'
user_col=f'new_subclass_{user}'
#save_dir = input("Please enter your working directory path: ")
save_path = save_path = f"{save_dir}/{user}_FRI_subclass_validation_test2.csv"

input_catalog='/beegfs/lofar/lclews/FRI_evolution/subclassification/inspection_results/new_env_final_classes_with_forces_and_pressures_model3.fits'


if os.path.exists(save_path):
    print("Previous classification file found,  resuming classification...")
    full_df =Table.read(save_path).to_pandas()
    
   
    # Ensure column exists AND purge NaNs
    if user_col not in full_df.columns:
        full_df[user_col] = ''
    else:
        full_df[user_col] = full_df[user_col].fillna('')
    
    # Enforce FITS-safe dtype
    full_df[user_col] = full_df[user_col].astype('<U15')
    
    # Now this works correctly
    df_to_do = df_to_do = full_df[
    (full_df['C_class'] == True) &
    (full_df['subclass_LC'] != 'unclassifiable') &
    (full_df['subclass_JC'] != 'unclassifiable') &
    (full_df['subclass_BM'] != 'unclassifiable') &
    (full_df[user_col] == '')]


  
    total_to_do = len(df_to_do)
    print(f"{len(df_to_do)} sources remaining for {user}.")

else:
    print("No previous file found  starting fresh.")
    print('Results will be saved to:', save_path)

    # Load original FITS catalogue
    base_table = Table.read(input_catalog)
    full_df = base_table.to_pandas()
    
    # Add new column for this user
         
    if user_col not in full_df.columns:
        full_df[user_col] = np.full(len(full_df), '', dtype='<U15')

   
    # Filter: C_class = True and no classifier said "unclassifiable"
    df_to_do = full_df[
    (full_df['C_class'] == True) &
    (full_df['subclass_LC'] != 'unclassifiable') &
    (full_df['subclass_JC'] != 'unclassifiable') &
    (full_df['subclass_BM'] != 'unclassifiable') &
    (full_df[user_col] == '')]

    
    
    
    if old_col in full_df.columns and user_col in full_df.columns:
        cols = list(full_df.columns)
    
        cols.remove(user_col)
        insert_at = cols.index(old_col) + 1
        cols.insert(insert_at, user_col)
    
        full_df = full_df[cols]


    # Save initial FITS
    Table.from_pandas(full_df).write(save_path, overwrite=True)
   
    print('starting visual inspection of %s sources' %len(df_to_do))
    total_to_do = len(df_to_do)
    
for index, row in df_to_do.iterrows():
    try:
        source = row['Source_Name'].decode() if isinstance(row['Source_Name'], bytes) else row['Source_Name']
        print('classifying source:', source)
        hostxdeg=row['optRA']
        hostydeg=row['optDec']
       
        description=classification(source, acceptable_response, hostxdeg,hostydeg)
     
        full_df.loc[index, user_col] = description
        full_df[user_col] = full_df[user_col].astype('<U15')
        Table.from_pandas(full_df).write(save_path, overwrite=True)
        
    except KeyboardInterrupt:
        print("\n Interrupted. Saving progress...")
        full_df[user_col] = full_df[user_col].astype('<U15')
        Table.from_pandas(full_df).write(save_path, overwrite=True)
        done = total_to_do - len(df_to_do.loc[full_df[user_col] == ''])
        total = total_to_do
        print(f"Progress: {done}/{total} sources completed "f"({100*done/total:.1f}%)")
        break
      

print('finished')




