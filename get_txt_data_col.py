# -*- coding: utf-8 -*-
"""
Created on Thu Oct 29 14:04:54 2020

@author: sabri
"""

""" 
This script defines 3 functions that will read in and sort SN data from .txt files, NOT FITS.
***Use this for getting data from individual files, use get_data to get whole epoch data sets. 
It was created based off the first SN2012au data recieved in the summer of 2020.
"""

from pathlib import Path
import re 
import numpy as np


def get_flux(folder_path, data_file, wave_start, flux_limit):
    """This function takes a .txt SN flux file with wavelength in column1 and flux in column2
    and returns only this data as [[wavelengths list], [flux lists]]. 
    The wave_start is the minimum wavelength value and the flux_limit is a conservative maximum flux value.
    Note: check that first and last few values of each list mathc the file 
    because values from the header can sneak in at the start and must be removed."""
    
    folder = Path(folder_path)
    file = folder/ data_file
    f = (open(file)).readlines() #read data line by line   
    column1=[] #wavelengths
    column2=[] #fluxs
    for i in f:
        lis = i.split()
        for n in lis:
            try:
                x = float(n)
                if x > wave_start: 
                    column1.append(x) #sorts wavlengths by setting lower threshold at wave_start
                elif 0 < np.abs(x) < flux_limit:
                    column2.append(x) #sort fluxes by seeting bounds between 0 and the flux_limit
                    
            except ValueError: #ignores the header words that cannot be made into floats
                pass
     
    if len(column1) != len(column2): #check to see if lengths match
        print("WARNING: data list lengths do not match: (" + str(len(column1)) + " vs " \
              + str(len(column2)) + ") Check lists against file.")
    return(column1, column2)

"""
data = get_flux("OG_data", "sn2012au_first2comb.flx.txt", 3000, .000001)
print(len(data[0]), len(data[1]))
#example of removing unwanted "wavelength" values from header
wavelength = data[0]
wavelength = wavelength[2:]
flux = data[1]
print(len(wavelength), len(flux))
print(wavelength[0:3], wavelength[998:])


data = get_flux("Epoch_1//all_comb", "sn2012au_allcomb.flx.txt", 3000, .000001)
print(len(data[0]), len(data[1]))
wavelength = data[0]; wavelength = wavelength[2:]
flux = data[1]
print(len(wavelength), len(flux))

"""
def get_qu(folder_path, data_file, wave_start, QU_min, QU_max):
    """This function takes a .txt SN q or u polarization file with wavelength in column1 and 
    q or u values in column2 and returns only this data as [[wavelengths list], [q or u lists]]. 
    The wave_start is the minimum wavelength value and the QU_min or max are conservative bounds 
    for these values. 
    Note: check that first and last few values of each list mathc the file 
    because values from the header can sneak in at the start and must be removed."""
    
    folder = Path(folder_path)
    file = folder/ data_file
    f = (open(file)).readlines() #read data line by line   
    column1=[] #wavelengths
    column2=[] #fluxs
    for i in f:
        lis = i.split()
        for n in lis:
            try:
                x = float(n)
                if x > wave_start: 
                    column1.append(x) #sorts wavlengths by setting lower threshold at wave_start
                elif QU_min < np.abs(x) < QU_max:
                    column2.append(x) #sort q or u by seeting bounds between given min and max
                    
            except ValueError: #ignores the header words that cannot be made into floats
                pass
     
    if len(column1) != len(column2): #check to see if lengths match
        print("WARNING: data list lengths do not match: (" + str(len(column1)) + " vs " \
              + str(len(column2)) + ") Check lists against file.")
    return(column1, column2)

"""
data = get_qu("OG_data", "sn2012au_first2comb.u.txt", 3000, -10, 10)
#example of removing unwanted "q" values from header
wavelength = data[0]
wavelength = wavelength[2:]
q = data[1]
q = q[21:]
print(len(wavelength), len(q))
print(wavelength[0:3], wavelength[998:])
print(q[0:25], q[998:])
"""

def get_qu_sig(folder_path, data_file, wave_start, error_max):
    """This function takes a .txt SN q or u polarization error files with wavelength in column1 
    and the error values in column2, and returns only this data as [[wavelengths list], [errors lists]]. 
    The wave_start is the minimum wavelength value and the error_max is a conservative maximum error value.
    Note: check that first and last few values of each list mathc the file 
    because values from the header can sneak in at the start and must be removed."""
    
    folder = Path(folder_path)
    file = folder/ data_file
    f = (open(file)).readlines() #read data line by line   
    column1=[] #wavelengths
    column2=[] #fluxs
    for i in f:
        lis = i.split()
        for n in lis:
            try:
                x = float(n)
                if x > wave_start: 
                    column1.append(x) #sorts wavlengths by setting lower threshold at wave_start
                elif 0 < np.abs(x) < error_max:
                    column2.append(x) #sort errors by seeting bounds between 0 and error_max
                    
            except ValueError: #ignores the header words that cannot be made into floats
                pass
     
    if len(column1) != len(column2): #check to see if lengths match
        print("WARNING: data list lengths do not match: (" + str(len(column1)) + " vs " \
              + str(len(column2)) + ") Check lists against file.")
    return(column1, column2)

"""    
data = get_qu_sig("OG_data", "sn2012au_first2comb.usig.txt", 3000, 1)
print(len(data[0]), len(data[1]))
wave = data[0]; wave = wave[2:]
print(len(wave))
"""