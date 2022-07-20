# -*- coding: utf-8 -*-
"""
Created on Sun Jan  3 14:15:41 2021

@author: sabri
"""
def Bin_data(data, old_binSize, binSize, z = 0.004483, new_filename = "none"):
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
    from get_data import P_tot, P_err
    import numpy as np
    from astropy.table import Table
    from astropy.io import ascii
    import os
    
    #access unbinned data from input table and de-redshift
    wave = data['wave'][:]/(1+ z)
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
    P = P_tot(qBinned, uBinned)
    Perr = P_err(P, qBinned, uBinned, qerrBinned, uerrBinned)
    
    binned = Table([waveBinned, fluxBinned, qBinned, qerrBinned, qsumBinned, uBinned, uerrBinned, P, Perr], names=['wave', 'flx', 'q', 'qerr', 'qsum', 'u', 'uerr', 'p', 'perr'])
    
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
    

