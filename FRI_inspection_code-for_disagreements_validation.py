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

def classification(
        source,
        acceptable_response,
        hostxdeg,
        hostydeg,
        validation_class,
        morph_class):
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
    plt.title(
    f"{source}\n"
    f"Validation: {validation_class}    Morph: {morph_class}",
    fontsize=10
    )
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
save_dir = input("Please enter your working directory path: ")
save_path = save_path = f"{save_dir}/{user}_FRI_subclass_validation_test2.csv"

input_catalog='/beegfs/lofar/lclews/FRI_evolution/subclassification/inspection_results/LC_FRI_subclass_validation.csv'


if os.path.exists(save_path):

    print("Previous classification file found, resuming classification...")
    full_df = pd.read_csv(save_path)

else:

    print("No previous file found, starting fresh.")
    full_df = pd.read_csv(input_catalog)

    if user_col not in full_df.columns:
        full_df[user_col] = ''

    full_df.to_csv(save_path, index=False)

# ----------------------------------------------------------
# get current Morph classifications from full catalogue
# ----------------------------------------------------------

full_catalog = Table.read(
    "/beegfs/lofar/lclews/FRI_evolution/subclassification/"
    "inspection_results/new_env_final_classes_with_forces_and_pressures_model3.fits"
).to_pandas()

full_catalog["Source_Name"] = (
    full_catalog["Source_Name"]
    .astype(str)
    .str.strip()
)

full_catalog["Morph"] = (
    full_catalog["Morph"]
    .astype(str)
    .str.strip()
)

full_df["Source_Name_1"] = (
    full_df["Source_Name_1"]
    .astype(str)
    .str.strip()
)

full_df["subclass"] = (
    full_df["subclass"]
    .astype(str)
    .str.strip()
)

comparison = full_df.merge(
    full_catalog[["Source_Name", "Morph"]],
    left_on="Source_Name_1",
    right_on="Source_Name",
    how="inner"
)

comparison = comparison[
    (comparison["subclass"].str.lower() != "unclassifiable") &
    (comparison["Morph"].str.lower() != "unclassifiable") &
    (comparison["subclass"] != comparison["Morph"])
].copy()

print(f"Found {len(comparison)} disagreement sources")

# make lookup dictionaries
validation_lookup = dict(
    zip(comparison["Source_Name_1"], comparison["subclass"])
)

morph_lookup = dict(
    zip(comparison["Source_Name_1"], comparison["Morph"])
)

sources_to_do = set(comparison["Source_Name_1"])

if user_col not in full_df.columns:
    full_df[user_col] = ''

full_df[user_col] = full_df[user_col].fillna('')

df_to_do = full_df[
    full_df["Source_Name_1"].isin(sources_to_do)
].copy()

df_to_do = df_to_do[
    df_to_do[user_col] == ''
]

total_to_do = len(df_to_do)

print(f"{total_to_do} sources remaining for {user}")
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
        validation_class = validation_lookup[source]
        morph_class = morph_lookup[source]
        
        description = classification(
            source,
            acceptable_response,
            hostxdeg,
            hostydeg,
            validation_class,
            morph_class
        )
             
        full_df.loc[index, user_col] = description
        full_df[user_col] = full_df[user_col].astype('<U15')
        full_df.to_csv(save_path, index=False)
        
    except KeyboardInterrupt:
        print("\n Interrupted. Saving progress...")
        full_df[user_col] = full_df[user_col].astype('<U15')
        full_df.to_csv(save_path, index=False)
        done = total_to_do - len(df_to_do.loc[full_df[user_col] == ''])
        total = total_to_do
        print(f"Progress: {done}/{total} sources completed "f"({100*done/total:.1f}%)")
        break
      

print('finished')




