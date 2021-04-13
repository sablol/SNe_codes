# -*- coding: utf-8 -*-
"""
Created on Sun Jan  3 14:15:41 2021

@author: sabri
"""

#importing tools
import pylab
import csv
import numpy as np
from astropy.table import Table
from astropy.io import ascii
import seaborn as sns
import scipy as sp
import astropy.io.fits as fits
import os
import glob
from astropy.io import fits as pyfits
import matplotlib.pyplot as plt
from get_data import *


def get_fits(file_path, file_pattern, new_filename = "none"):
    """Takes in the filepath to a folder that contains .fits files
    NOTE: make sure the files in the folder are in the order specified by names:
     "wave", "q", "qerr", "qsum", "u", "uerr", "usum"  
    as that is the order they will be written into the table. 
    see astropy.table accessing table to read data out of table""" 

    #file_pattern = "\*.fits" 
    fileList = glob.glob(file_path + file_pattern)
   
    data_array = [] 

    for file in fileList: 
        data = pyfits.open(file)
        temp_list = []
        
        for i in data[0].data:
            temp_list.append(i) #append individual values to temp list
        data_array.append(temp_list) #append temp list to data array 
    
    #Calculate wavelength values based on header info
    length = float(data[0].header['NAXIS1'])
    start = float(data[0].header['CRVAL1'])
    step = float(data[0].header['CDELT1'])
    stop = start + (length*step)
    waves = np.arange(start,stop, step)
    data_array = [waves] + data_array #add list of waves to begining of data array
    
    output = Table(data_array, names = ["wave", "q", "qerr", "qsum", "u", "uerr", "usum"])
    if new_filename != "none":
        #optional if a new_filename is given write an output txt file
        output.write(new_filename+".txt", format = 'ascii')

    return(output)


def Bin_data(data, old_binSize, binSize, z, new_filename = "none"):
    """Function bins data and takes in five arguments. 
    The first is a data table with columns named: wave, qsum, q, qerr, usum, u, uerr
    *NOTE: columns must be in this order with these names in order to work.
    **output from get_fits can be used as data input
    old_binSize and binSize are single numbers.
    z is the redift for the SN or host galaxy.
    A new_filename is optional to save the ouput to, if blank no file is created. 
    Returns a table of binned data in the same format as the input data table. 
    """
    
    #access unbinned data from input table and de-redshift
    wave = data['wave'][:]/(1+ z)
    #print("wave type from input = ", type(wave))
    q = data['q']
    qerr = data['qerr']
    qsum = data['qsum'] 
    u = data['u']
    uerr = data['uerr']
    usum = data['usum']
    
    #combine flux with q/u 
    Q_flux = q*qsum
    Q_err_sq = (qerr*qsum)**2
    U_flux = u*usum
    U_err_sq =(uerr*usum)**2

    #bin wavelengths
    wmin = wave.min() 
    wmax = wave.max()
    #binSize = 10
    #old_binSize = 4
    bins = np.arange(wmin, wmax, binSize) #bin ranges based on data 
    inBin = np.digitize(wave, bins=bins, right = 'True') #what bin data falls in 
    binArray = np.arange(np.max(inBin)) #bins numbered 
    binCheck = (binArray[:, None] == inBin[None, :]) #True/False array: whether or not data falls in bin 
    waveBinned = (wave[None,:]*binCheck).sum(axis=1)/binCheck.sum(axis=1)
    #print("wave type from binning = ", type(waveBinned))
    #Check bin wavelengths is correct #print(len(wave), len(waveBinned)) #print(wave, inBin, waveBinned)

    #bin flux
    qfluxBinned = (qsum[None,:]*binCheck).sum(axis=1) #adding all flux values in each bin (relative so no need for average)
    ufluxBinned = (usum[None,:]*binCheck).sum(axis=1)

    #bin stokes
    qBinned = ((Q_flux[None,:]*binCheck).sum(axis=1))/(qfluxBinned) #q*qflux binned/qflux binned= q binned
    uBinned = ((U_flux[None,:]*binCheck).sum(axis=1))/(ufluxBinned) #same for u

    #bin stokes Errors
    binChange = np.sqrt(old_binSize/binSize)
    qerr_quad = np.sqrt((Q_err_sq[None,:]*binCheck).sum(axis=1)) #add qerr*qflux in quadrature
    qerrBinned = (qerr_quad*binChange)/qfluxBinned #factor in change in bin size and remove qflux 
    uerr_quad = np.sqrt((U_err_sq[None,:]*binCheck).sum(axis=1))#same for uerr
    uerrBinned = (uerr_quad*binChange)/ufluxBinned

    binned = Table([waveBinned, qfluxBinned, qBinned, qerrBinned, ufluxBinned, uBinned, uerrBinned], names=['wave', 'qsum', 'q', 'qerr', 'usum', 'u', 'uerr'])
    
    if new_filename != "none":
    #optional if a new_filename is given write an output txt file
        binned.write(new_filename+".txt", format = 'ascii')
        print("Data binned to " + str(binSize) + " Angstroms and saved to " + str(new_filename) +".txt")
    
    print("Data binned to " + str(binSize) + " Angstroms. To save output provide filename.")
    return(binned)
    
#test = get_txtFITS("Epoch_1\\all_comb", "sn2012au_allcomb.qsum.txt", "sn2012au_allcomb.q.txt", "sn2012au_allcomb.qsig.txt", "sn2012au_allcomb.qsum.txt", "sn2012au_allcomb.u.txt", "sn2012au_allcomb.usig.txt", "sn2012au_allcomb.unbinned.txt")
#print(test)
#bin_test = Bin_data(test, 4, 8, 0)    
#print(bin_test)
