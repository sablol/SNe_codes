# -*- coding: utf-8 -*-
"""
Created on Sun Jan  3 14:15:41 2021

@author: sabri
"""
def get_fits(file_path, names = ["wave", "q", "qerr", "qsum", "u", "uerr", "usum"], file_pattern = "\*.fits", new_filename = "none", z = 0):
    """Takes in the filepath to a folder that contains .fits files
    file_pattern = "\*.fits unless specified otherwise.
    NOTE: make sure the files in the folder are in the order specified by names list, the default is:
     "q", "qerr", "qsum", "u", "uerr", "usum"  
    as that is the order they will be written into the table with a calculated wave column first. 
    see astropy.table accessing table to read data out of table
    To use output with Bin_data usum column must be renamed flux.""" 

    import glob
    from astropy.io import fits as pyfits
    import astropy.io.fits as fits
    from astropy.table import Table
    from astropy.io import ascii
    import numpy as np

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
    length = float(data[0].header['NAXIS1']); 
    start = float(data[0].header['CRVAL1']); 
    step = float(data[0].header['CDELT1']); 
    #print(data[0].header['CDELT1'])
    stop = start + (length*step);
    waves = np.arange(start, stop, step);
    waves = waves[:]/(1+ z)
    #print(len(data_array[0]))

    if len(waves) == len(data_array[1]):
        print('waves length: ' + str(len(waves)) + ', data array column lengths: ' +str(len(data_array[1])))
        data_array = [waves] + data_array #add list of waves to begining of data array
        output = Table(data_array, names = names)
        if new_filename != "none":
            #optional if a new_filename is given write an output txt file
            output.write(new_filename+".txt", format = 'ascii')
        return(output)    
    
    if len(waves)-1 == len(data_array[1]):
        print('Warning: uneven list values, list length of wavelengths calculated from header is one longer than data array length.\n See get_fits function to review adjustment made for lists to be equal.')
        print('Originally calculated lengths: ')
        print('waves length: ' + str(len(waves)) + ', data array column lengths: ' +str(len(data_array[1])))
        waves = np.arange(start, stop-(step-1), step); 
        print('Recalculated lengths:')
        print('waves length: ' + str(len(waves)) + ', data array column lengths: ' +str(len(data_array[1]))) #if len(wave) is one longer than len(data) use this line instead
        data_array = [waves] + data_array #add list of waves to begining of data array
        output = Table(data_array, names = names)
        if new_filename != "none":
            #optional if a new_filename is given write an output txt file
            output.write(new_filename+".txt", format = 'ascii')
        return(output)
    
    else: 
        print('uneven list values, data array returned with [0] = wavelengths calculated from header')
        return(data_array)

def Bin_data(data, old_binSize, binSize, z = 0, new_filename = "none"):
    """Function bins data and takes in five arguments. 
    The first is a data table with columns named: wave, flx, q, qerr, qsum, u, uerr
    *NOTE: columns must be in this order with these names in order to work.
    **output from get_fits can be used as data input
    old_binSize and binSize are single numbers.
    z is the redift for the SN or host galaxy.
    A new_filename is optional to save the ouput to, if blank no file is created. 
    Returns a table of binned data in the same format as the input data table with added total polarization columns. 
    
    Example:
    data1 = get_txtFITS("Epoch_1/all_comb", "sn2012au_allcomb.qsum.txt", "sn2012au_allcomb.q.txt", "sn2012au_allcomb.qsig.txt", "sn2012au_allcomb.flx.txt", "sn2012au_allcomb.u.txt", "sn2012au_allcomb.usig.txt", "sn2012au_allcomb.unbinned.txt")
    epoch1 = Bin_data(data1, 4, 20, 0.004483, "text_doc_name.txt") 
    """
    #importing tools/functions
    from get_data import P_tot, P_err, P_debi3
    import numpy as np
    from astropy.table import Table
    from astropy.io import ascii
    import os
    
    #access unbinned data from input table and de-redshift
    wave = data['wave']#[:]/(1+ z) #uncomment if unbinned data is not deredshifted 
    #print("wave type from input = ", type(wave))
    q = data['q']*100
    qerr = data['qerr']*100
    qsum = data['qsum'] 
    u = data['u']*100
    uerr = data['uerr']*100
    flx = data['flx']
    
    
    #combine flux with q/u 
    Q_flux = q*flx
    Q_err_sq = (qerr*flx)**2
    U_flux = u*flx
    U_err_sq =(uerr*flx)**2

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
    fluxBinned = (flx[None,:]*binCheck).sum(axis=1) #adding all flux values in each bin (relative so no need for average)
    
    #bin qsum 
    qsumBinned = (qsum[None,:]*binCheck).sum(axis=1) #adding all qsum in each bin (I think they are relative too so not averaged)
    
    #bin stokes
    qBinned = ((Q_flux[None,:]*binCheck).sum(axis=1))/(fluxBinned) #q*qflux binned/qflux binned= q binned
    uBinned = ((U_flux[None,:]*binCheck).sum(axis=1))/(fluxBinned) #same for u

    #bin stokes Errors
    binChange = np.sqrt(old_binSize/binSize)
    qerr_quad = np.sqrt((Q_err_sq[None,:]*binCheck).sum(axis=1)) #add qerr*qflux in quadrature
    qerrBinned = (qerr_quad*binChange)/fluxBinned #factor in change in bin size and remove qflux 
    uerr_quad = np.sqrt((U_err_sq[None,:]*binCheck).sum(axis=1))#same for uerr
    uerrBinned = (uerr_quad*binChange)/fluxBinned
    
    #calculate total polarization and errors 
    P_OG = P_tot(qBinned, uBinned)
    Perr = P_err(P_OG, qBinned, uBinned, qerrBinned, uerrBinned)
    #P = P_debi(qBinned, uBinned, qerrBinned, uerrBinned)
    P = P_debi3(P_OG, Perr)
    
    binned = Table([waveBinned, fluxBinned, qBinned, qerrBinned, qsumBinned, uBinned, uerrBinned, P, Perr], names=['wave', 'flx', 'q', 'qerr', 'qsum', 'u', 'uerr', 'p', 'perr'])
    
    if new_filename != "none":
    #optional if a new_filename is given write an output txt file
        binned.write(new_filename+".txt", format = 'ascii', overwrite = True)
        #with open(new_filename+".txt", 'r+') as f:
            #content = f.read()
            #f.seek(0, 0)
            #f.write("This data has been deredshifted using a value of: "+ str(z))
        print("Data binned to " + str(binSize) + " Angstroms and saved to " + str(new_filename) +".txt")
    else:
        print("Data binned to " + str(binSize) + " Angstroms. To save output provide filename.")
    return(binned)

def BinFLX_CSV(file, binSize, z = 0, new_filename = 'none'):
    """Function bins non-SPOL spectra (only wavelength and flx). 
    file parameter should point directly to a single CSV file.
    binSiz should be specified in whole integers (old binsize is not necessary because errors on flx cannot be calculated)
    If a z value is given binned data will be deredshifted accordingly. 
    If a new_filename is provided a txt file of binned dat will be saved."""
    
    import pandas as pd
    import numpy as np
    from astropy.table import Table
    from astropy.io import ascii
    import os
    
    #open non-spol CSV data
    #for OG file from Stan
    #data = pd.read_csv(file, sep = ',')
    #wave = data.Wavelength.tolist() #for OG file from Stan
    #flx = data.Flux.tolist()#for OG file from Stan
    #data = Table([wave, flx], names = ['wave', 'flx']) #
    #wave = data['wave'][:]/(1+ z)
    #flx = data['flx']
    
    #for data converted to table given in function arguments  
    wave = file['wave'][:]/(1+ z)
    flx = file['flx']
    #bin wavelengths
    wmin = wave.min() 
    wmax = wave.max()
    bins = np.arange(wmin, wmax, binSize) #bin ranges based on data 
    inBin = np.digitize(wave, bins=bins, right = 'True') #what bin data falls in 
    binArray = np.arange(np.max(inBin)) #bins numbered 
    binCheck = (binArray[:, None] == inBin[None, :]) #True/False array: whether or not data falls in bin 
    waveBinned = (wave[None,:]*binCheck).sum(axis=1)/binCheck.sum(axis=1)
    #bin flux
    fluxBinned = (flx[None,:]*binCheck).sum(axis=1) #adding all flux values in each bin (relative so no need for average)
    binned = Table([waveBinned, fluxBinned], names=['wave', 'flx'])

    #save file 
    if new_filename != "none":
    #optional if a new_filename is given write an output txt file
        binned.write(new_filename+".txt", format = 'ascii')
        #with open(new_filename+".txt", 'r+') as f:
            #content = f.read()
            #f.seek(0, 0)
            #f.write("This data has been deredshifted using a value of: "+ str(z))
        print("Data binned to " + str(binSize) + " Angstroms and saved to " + str(new_filename) +".txt")
    else:
        print("Data binned to " + str(binSize) + " Angstroms. To save output provide filename.")
    return(binned)    

