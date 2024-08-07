# -*- coding: utf-8 -*-
"""
Created on Mon Jul 12 10:41:49 2021

@author: sabri
"""
#import libraries used
import numpy as np
import spectres as spec
from astropy.table import Table
from get_data import *

def spec_Bin_data(data, binSize, z, new_filename = "none"):
    """Function bins data and takes in four arguments. 
    The first is a data table with columns named: wave, flx, q, qerr, qsum, u, uerr
    *NOTE: columns must be in this order with these names in order to work.
    **output from get_txtFITS or get_fits (after renaming column) can be used as data input
    
    binSize is how many Angstroms to include for each data point
    z is the redift for the SN or host galaxy.
    A new_filename is optional to save the ouput to, if blank no file is created. 
    Returns a table of binned data in the same format as the input data table with added total polarization columns.
    Returns Q & U data as percentages!!!!
    
    Examples:
    data1 = get_txtFITS("Epoch_1/all_comb", "sn2012au_allcomb.qsum.txt", "sn2012au_allcomb.q.txt", "sn2012au_allcomb.qsig.txt", "sn2012au_allcomb.flx.txt", "sn2012au_allcomb.u.txt", "sn2012au_allcomb.usig.txt")
    epoch1 = Bin_data(data1, 20, 0.004483, "new_binned_datafile_name") 
    data6 = get_fits("Epoch_6"); data6['usum'].name = 'flx'
    epoch6 = Bin_data(data6, 20, 0.004483)
    """
    
    #access unbinned data from input table and de-redshift
    wave = data['wave'][:]/(1+ z)
    q = data['q']*100
    qerr = data['qerr']*100
    qsum = data['qsum'] 
    u = data['u']*100
    uerr = data['uerr']*100
    flx = data['flx']

    #specify binning edges
    wave_bins = np.arange(wave[0], wave[-1], binSize)

    #bin flux 
    flux_binned = spec.spectres(wave_bins, wave, flx)

    #bin Q/U (multiply by flux first, then bin, then divide flux out)
    q_binned_data = (spec.spectres(wave_bins, wave, q*flx, spec_errs=(qerr*flx)))/flux_binned
    q_binned = q_binned_data[0]; qerr_binned = q_binned_data[1]
    qsum_binned_data = spec.spectres(wave_bins, wave, qsum, spec_errs=(qerr)) #binned like a flux not a polarization vector 
    qsum_binned = qsum_binned_data[0]
    u_binned_data = (spec.spectres(wave_bins, wave, u*flx, spec_errs=(uerr*flx)))/flux_binned
    u_binned = u_binned_data[0]; uerr_binned = u_binned_data[1]

    #calculate total polarization and errors (p_tot function in get_data)
    P = P_tot(q_binned, u_binned)
    Perr = P_err(P, q_binned, u_binned, qerr_binned, uerr_binned)
    
    #create table of binned data to return
    binned_data = Table([wave_bins, flux_binned, q_binned, qerr_binned, qsum_binned, u_binned, uerr_binned, P, Perr], names=['wave', 'flx', 'q', 'qerr', 'qsum', 'u', 'uerr', 'p', 'perr'])

    #optional if a new_filename is given write an output txt file    
    if new_filename != "none":
        binned_data.write(new_filename+".txt", format = 'ascii')
        print("Data binned to " + str(binSize) + " Angstroms and saved to " + str(new_filename) +".txt")
    else:
        print("Data binned to " + str(binSize) + " Angstroms. To save output provide filename.")
    
    return(binned_data)
