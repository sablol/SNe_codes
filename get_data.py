# -*- coding: utf-8 -*-
"""
Created on Mon Nov  9 17:14:01 2020

@author: sabri
"""

"""
Simple way to read in binned data: 
from astropy.io import ascii
from astropy.table import table
data = astropy.io.ascii.read(file, delimiter='\s', guess = False, names =  ['waves','qsum', 'q', 'qerr', 'usum', 'u', 'uerr'])
print(type(data['waves'][5]))
""" 
import numpy as np
from scipy.odr import *
import scipy.stats as stats

def z_correct(wave, z):
    """takes in wavelength (or list of wavelengths) and redshift value, corrects for it and returns new wavelength(s)"""  
    wave = wave/(1+ z)
    return(wave)

def get_binned(filename, z = 0):
    """Takes in binned data folder/filename and returns a table of data with wavelengths corrected for redshift if redshift is provided. 
       Make sure input file is in the order of wave, qsum, q, qerr, usum, u, uerr and the return table will be ordered the same.
       Example: data = get_binned("Epoch_1/old_binning/sn2012au_allcomb_b10.txt", 0.004483)""" 
    
    from astropy.io import ascii
    data = ascii.read(filename, delimiter='\s', guess = False, names =  ['wave','qsum', 'q', 'qerr', 'usum', 'u', 'uerr'])
    
    if z != 0:
        data['wave'][:] = z_correct(data['wave'][:], z)
    
    return(data)   

def P_tot(Q, U):
    """Calculates total Polarization in quadrature from Q and U values.
    Data is returned as an astropy table column."""
    from astropy.table import Table
    import numpy as np
    
    P = (np.sqrt(np.power(Q, 2) + np.power(U, 2)))
    t = Table([P], names = ['p'])
    return(t['p'])

def P_err(P, Q, U, Qerr, Uerr):
    """Calculates total Polarization errors in quadrature from Q, Q errors and U, U errors values.
    Data is returned as an astropy table column."""
    from astropy.table import Table
    import numpy as np
    #Wardle Kronberg 1974 sig_p ~ sig_p_debi
    #used in 2010jl paper
    Perr = (1/P)*(np.sqrt(np.power(Q, 2)*np.power(Qerr, 2) + np.power(U, 2)*np.power(Uerr, 2)))
    t = Table([Perr], names = ['perr'])
    return(t['perr'])

def P_debi(Q, U, Qerr, Uerr):
    from astropy.table import Table
    import numpy as np 
    #Christopher Bilinski 2018 ref. Wardle Kronberg 1974 
    P = (np.sqrt(np.power(Q, 2) + np.power(U, 2) - 0.5*(np.power(Qerr, 2) + np.power(Uerr, 2))))
    t = Table([P], names = ['p'])
    return(t['p'])

def P_debi2(P, Perr):
    from astropy.table import Table
    import numpy as np
    #Wardle Kronberg 1974 (2010jl for cases where P_sig > p)
    Pdebi = []
    for i in range(len(P)):
        radicand = P[i]**2 - Perr[i]**2
        if radicand > 0 :
            Pdebi.append(+(np.sqrt(np.absolute(radicand))))
        else:
            Pdebi.append(-(np.sqrt(np.absolute(radicand))))
    t = Table([Pdebi], names = ['p'])
    return(t['p'])

def P_debi3(P, Perr):
    from astropy.table import Table
    #Wardle Kronberg 1974 used in 2010jl paper
    #used for 2012au paper
    Pdebi = []
    for i in range(len(P)):
        if Perr[i] > P[i]: #use P_debi2 in this case 
            radicand = P[i]**2 - Perr[i]**2
            if radicand > 0 :
                Pdebi.append(+(np.sqrt(np.absolute(radicand))))
            else:
                Pdebi.append(-(np.sqrt(np.absolute(radicand))))
        else:
            Pdebi.append(P[i]*np.sqrt(1-(Perr[i]/P[i])**2))
    t = Table([Pdebi], names = ['p'])
    return(t['p'])

def QRSP(Q, U, Theta):
    """Calculates the Rotated Stokes Parameter (RSP) version of P given a list of Q and U values and a single theta value. 
    (Theta should be smoothed or relatively constant with wavelength)
    Data is retruned as a astropy table column."""
    from astropy.table import Table
    import numpy as np
    QRSP = (Q*np.cos(2*np.radians(Theta)) + U*np.sin(2*np.radians(Theta)))
    t = Table([QRSP], names = ['qrsp'])
    return(t['qrsp'])

def URSP(Q, U, Theta):
    """Calculates the Rotated Stokes Parameter (RSP) version of P given a list of Q and U values and a single theta value. 
    (Theta should be smoothed or relatively constant with wavelength)
    Data is retruned as a astropy table column."""
    from astropy.table import Table
    import numpy as np
    URSP = (-Q*np.sin(2*np.radians(Theta)) + U*np.cos(2*np.radians(Theta)))
    t = Table([URSP], names = ['ursp'])
    return(t['ursp'])

def get_smooth_rsp(data, smooth_window, rsp_stokes = 'qrsp'):
    """Calculates the Rotated Stokes Parameter (RSP) version of P given a data table containing q and u columns 
    and smoothing window (float value in Å) to calculate the average PA within to use for smoothing angle in QRSP or URSP function. 
    Data table given as input is retruned with an added table column data['qrsp']*
    *Note: ursp is calculated and returned instead if rsp_stokes is specified as something other than 'qrsp'. """
    
    sm_wave_starts = np.arange(start = data['wave'][0], stop = data['wave'][-1], step = smooth_window)
    sm_theta = []
    rsp = []
    for wave_start in sm_wave_starts:
        sm_reg = get_lines(data, wave_start, wave_start + smooth_window)
        ave_q = get_weighted_mean(sm_reg)['q']
        ave_u = get_weighted_mean(sm_reg)['u']
        ave_theta = PA(ave_q, ave_u); sm_theta.append(ave_theta)
        if rsp_stokes == 'qrsp':
            rsp = np.hstack([rsp, QRSP(sm_reg['q'], sm_reg['u'], ave_theta)])
        else:
            rsp = np.hstack([rsp, URSP(sm_reg['q'], sm_reg['u'], ave_theta)])
    data[rsp_stokes]= rsp
    return(data, sm_wave_starts, sm_theta)


def PA(Q, U, adjust_pa = 'off'):   
    from astropy.table import Table
    import numpy as np
    """Calculates PA and adjusts calculated value by adding 180 to all PA's if PA = 'on' 
    or adds 180 to PA's < adjust_pa value specified. Default is adjust_pa = 'off' so PA returned
    is the raw calculated value. Returns a table with column PA. Used for get_PA_column."""
    PA = np.degrees(.5*np.arctan2(U, Q)) #calculate PA
    
    for i in range(len(PA)):
        if adjust_pa == 'on':  #add 180 to all PA's if adjust = 'on' 
            PA[i] = PA[i]+180
        elif adjust_pa == 'off': #leave all PA as is 
            PA[i] = PA[i] 
        elif adjust_pa == 'compare':
           while abs(PA[i] - PA[i-1]) > 90: #add 180 if the difference between two consecutive PA's is greater than 90
               PA[i] = PA[i] + 180
        else: 
            if PA[i] < adjust_pa: #add 180 to PA's less than specified adjust_pa value
                PA[i] = PA[i]+180

    t = Table([PA], names = ['PA'])
    return(t['PA'])

def PA_err(P, Q, U, Qerr, Uerr):
    """Calculates total Polarization Angle errors in quadrature from Q, Q errors and U, U errors values and returns errors as degrees.
    Data is returned as an astropy table column."""
    from astropy.table import Table
    import numpy as np
    
    PAerr = (1/P**2)*np.sqrt((Q**2)*(Uerr**2) + (U**2)*(Qerr**2))*(90/np.pi)
    #PAerr = (1/P**2)*np.sqrt((qspec_resample**2)*(uspec_errs_resample**2) + (uspec_resample**2)*(qspec_errs_resample**2))*(180/np.pi)   
    t = Table([PAerr], names = ['pa_err'])
    return(t['pa_err'])

def get_pa_column(data, adjust = 'off'):
    """Takes in data table and adds a pa and pa error column. 
    Does not adjust PA unless adjust = 'on' is specified. """
    PAs = PA(data['q'], data['u'], adjust_pa = adjust)
    PA_errs = PA_err(data['p'], data['q'], data['u'], data['qerr'], data['uerr'])

    data['pa'] = PAs
    data['paerr'] = PA_errs
    return(data)

def get_QU(P, PA):
    import numpy as np
    q = P*np.cos(2*np.radians(PA))
    u = P*np.sin(2*np.radians(PA))
    return(q, u)


print(get_QU([0.22], [53.03]))

def skip_to(file, line, col2_name, **kwargs):
    import os
    from astropy.table import Table
    
    """Function purpose: read in a txt file formated with a fits header and return data without the header
    file is the .txt file 
    line is the last line before the data starts (usually "END" for fits) 
        *NOTE: for skip_to to work there has to be one blank line between the line specified and the data, no more, no less. 
        If code doesn't work check the file provided for the spaced between header and data. 
    col2_name is what you would like to call the second column of data in the file given
    **kwargs specifies the sepperation delimiter usually sep="\s+"
    Example: #data_no_header = skip_to("C:\\Users\\sabri\\Desktop\\SN2012_au\\Epoch_1\\sn2012au_allcomb.flx.txt","END", "qsum", sep="\s+")
    """
    
    wave = []
    col2 = []
    if os.stat(file).st_size == 0:
        raise ValueError("File is empty")
    with open(file) as f:
        pos = 0
        cur_line = f.readline()
        while not cur_line.startswith(line):
            pos = f.tell()
            cur_line = f.readline()
        f.seek(pos)
        
        f.readline() #skips line with "END"
        f.readline() #skips blank line after "END" line - assumes only one blank line after 
        
        for stuff in f:
            
            item = stuff.split()
            wave.append(float(item[0]))
            col2.append(float(item[1]))
            
        return Table([wave, col2], names = ['wave', str(col2_name)])

def get_txt_data(dir, file_string, new_filename):
    """Gets data from two column format txt files where first column is wavelengths and the next is 
    flx, q, qsig, qsum, u, usig
    Files are called by the directory specified then the file_string that is the common name of each file and then adds the .flx.txt
    extension. *This will not work if they are not saved with the second column of data followed by the .txt extension.
    It then combines data into a new data table including polarization and polarization error values saved to the new_filename 
    specified and in the same folder as old data. Output of this can then be used with Bin_Data() function"""
    
    import pylab
    from astropy.table import Table

    wave = pylab.loadtxt(dir+file_string+'flx.txt')[:,0] 
    flx = pylab.loadtxt(dir+file_string+'flx.txt')[:,1]
    qsum = pylab.loadtxt(dir+file_string+'qsum.txt')[:,1]
    q = pylab.loadtxt(dir+file_string+'q.txt')[:,1]
    qerr = pylab.loadtxt(dir+file_string+'qsig.txt')[:,1]
    u = pylab.loadtxt(dir+file_string+'u.txt')[:,1]
    uerr = pylab.loadtxt(dir+file_string+'usig.txt')[:,1]
    
    #compute p values
    p = P_tot(q, u)
    perr = P_err(p, q, u, qerr, uerr)
    
    output_data = Table([wave, flx, q, qerr, qsum, u, uerr, p, perr], names = ['wave', 'flx', 'q', 'qerr', 'qsum', 'u', 'uerr', 'p', 'perr'])
    if new_filename != "none":
    #optional if a new_filename is given write an output txt file
        output_data.write(dir+ new_filename+".txt", format = 'ascii', overwrite = True)
    
    return(output_data)

def get_txtFITS(folder_path, flxfile, qfile, qerrfile, qsumfile, ufile, uerrfile, z = 0.004483, new_filename = "none"):    
    """
    **NOTE: this works for txt files made exactly from FITS files and uses the fits 
    headings so if the file is not the txt version of the fits this will not work!**
    Takes in sepperate txt files for qsum, q, qerr, usum, u, uerr data (wavelength is the first column in all of them)
    and returns one table with all of it combined, if new_filename is 
    provided this table is saved in a .txt file in the epoch folder
    test = get_txtFITS("Epoch_1/all_comb", "sn2012au_allcomb.qsum.txt", "sn2012au_allcomb.q.txt", "sn2012au_allcomb.qsig.txt", "sn2012au_allcomb.qsum.txt", "sn2012au_allcomb.u.txt", "sn2012au_allcomb.usig.txt", "sn2012au_allcomb.unbinned.txt")
    Table output can then be fed to Bin_Data() function for binning. 
    """
    
    from pathlib import Path
    import re 
    import numpy as np
    from astropy.table import Table
    from astropy.io import ascii
    
    #The first is a data table with columns named: wave, flx, q, qerr, qsum, u, uerr
    
    folder = Path(folder_path)    
#get wavelengths and flux values
    file1 = folder/ flxfile
    data1 = skip_to(file1, "END", "flx", sep="\s+")
    wave = data1['wave'][:]/(1+ z)
    flx = data1['flx']
#q values    
    file2 = folder/ qfile
    data2 = skip_to(file2, "END", "q", sep="\s+")
    q = data2['q']
#qerr values    
    file3 = folder/ qerrfile
    data3 = skip_to(file3, "END", "qerr", sep="\s+")
    qerr = data3['qerr']
#usum values     
    file4 = folder/ qsumfile
    data4 = skip_to(file4, "END", "qsum", sep="\s+")
    qsum = data4['qsum']
    
#u values     
    file5 = folder/ ufile
    data5 = skip_to(file5, "END", "u", sep="\s+")
    u = data5['u']

#uerr values     
    file6 = folder/ uerrfile
    data6 = skip_to(file6, "END", "uerr", sep="\s+")
    uerr = data6['uerr']  
    
#compute p values
    p = P_tot(q, u)
    perr = P_err(p, q, u, qerr, uerr)
    
    output_data = Table([wave, flx, q, qerr, qsum, u, uerr, p, perr], names = ['wave', 'flx', 'q', 'qerr', 'qsum', 'u', 'uerr', 'p', 'perr'])
    if new_filename != "none":
    #optional if a new_filename is given write an output txt file
        output_data.write(str(folder)+"/"+ new_filename+".txt", format = 'ascii')
    
    return(output_data)
 

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
                 
    
def get_lines(data, minn, maxx):
    from astropy.table import Table
    from plot_data import P_tot
    import numpy as np
    """gets data from lines regions from data tables with columns labeled as:
    wave, flx, q, qerr, u, uerr, p, perr
    based on min and max wavelengths provided. Returns a table with columns named the same.
    It will work for whatever columns are present, the table does not need to have all of them."""
    
    """more efficient way:
        indeces = np.where((epoch2_20['wave'] > 6000) & (epoch2_20['wave'] < 7000))
        print(epoch2_20[indeces])"""
    
    wl = []
    flx = []
    q = []
    qerr = []
    u = []
    uerr = []
    p = []
    perr = []
    
    data = data.to_pandas()
    line_data = Table()
    
    for i in range(0, len(data["wave"])):
        w = float(data["wave"][i])
        if w >= minn and w <= maxx:
            wl.append(data['wave'][i])
            if 'flx' in data:
                flx.append(data['flx'][i])
            if 'q' in data: 
                q.append(data['q'][i])
            if 'qerr' in data: 
                qerr.append(data['qerr'][i])
            if 'u' in data: 
                u.append(data['u'][i])
            if 'uerr' in data: 
                uerr.append(data['uerr'][i])
            if 'p' in data: 
                p.append(data['p'][i])
            if 'perr' in data: 
                perr.append(data['perr'][i])
            
    line_data['wave'] = wl
    if len(flx) == len(wl):
        line_data['flx'] = flx
    if len(q) == len(wl):
        line_data['q'] = q
    if len(qerr) == len(wl):
        line_data['qerr'] = qerr
    if len(u) == len(wl):
        line_data['u'] = u
    if len(uerr) == len(wl):
        line_data['uerr'] = uerr
    if len(p) == len(wl):
        line_data['p'] = p
    if len(perr) == len(wl):
        line_data['perr'] = perr
    return(line_data)   



def get_redshift_velocity(z):
    c = 299792
    v = c*z
    return(v)

def get_velocity(rest, obs):
    """Calculates velocity from rest wavelength and observed wavelength, returns velocity in km/s. 
    Takes in two inputs, the restwavelength and the observed wavelength. 
    """
    conversion = 10**-12 #angstroms to km
    c = 299792
    
    return(((obs*conversion-rest*conversion)/(rest*conversion))*c)

def get_vel_column(data, rest_wave):
    """Takes in data table and rest wavelength for line of interest and calculates velocity 
    with respect to that wavelength. Useful with output of get_lines and to plot velocity space."""
    vel_list =[]
    for i in data['wave']:
        vel_list.append(get_velocity(rest_wave, i))
    
    data['vel'] = vel_list
    return(data)

def get_obswave(rest, velocity):
    """Calculates the observed wavelength both red and blue shifted from a given rest wavelength 
    and velocity. Returns (observed blue wavelength, observed red wavelength) """
    c = 299792
    red_obs = rest - (velocity/c)*rest
    blue_obs = rest + (velocity/c)*rest
    return(blue_obs, red_obs)

def get_aves(epoch_list, region):
    """Function takes a list of epoch data tables (from Bin_data output) and a list of two wavelengths
    that it will then calculated the average value of everything in between using get_lines. 
    Returned is a data table with column labels: ['wave', 'flx', 'q', 'qerr', 'u', 'uerr', 'p', 'perr']
    the length of the table is the number of epochs in the epoch_list provided and each value is the 
    mean or error-weighted mean for that data type and epoch.   
    ex. input: cont1 = get_aves(epoch_data_table, [6400, 6700])
    """
    import numpy as np
    from astropy.table import Table
    from statsmodels.stats.weightstats import DescrStatsW
           
    wave = []; flx = []; q = []; qerr = []; u =[]; uerr =[]; p =[]; perr =[] #lists to store averages 

    for epoch in epoch_list:
        lines = get_lines(epoch, region[0], region[1]) #get only data for desired regions 
        #calculate averages **Note for list x, x[~numpy.isnan(x)] gets rid of nans in list and is used for each data list
        wave_ave = np.nanmean(lines['wave'][~np.isnan(lines['wave'])]); wave.append(wave_ave)
        flx_ave = np.nanmean(lines['flx'][~np.isnan(lines['flx'])]); flx.append(flx_ave)
        
        #for polarization calculate error-weighted mean and error on the wieghted mean 
        #Note: there is a list that holds all mean values and a list for q, u, qerr and uerr that updates with  
        #a single value each time b/c the function to calculate P and Perr take lists as inputs
        q_stats = DescrStatsW(lines['q'][~np.isnan(lines['q'])], weights=lines['qerr'][~np.isnan(lines['qerr'])], ddof=0)
        q.append(q_stats.mean); q_mean = [q_stats.mean] 
        q_mean_err = (np.sqrt(1/np.sum(1/(lines['qerr'][~np.isnan(lines['qerr'])])**2)))
        qerr.append(q_mean_err); #q_mean_errlist = [q_mean_err]
        
        u_stats = DescrStatsW(lines['u'][~np.isnan(lines['u'])], weights=lines['uerr'][~np.isnan(lines['uerr'])], ddof=0)
        u.append(u_stats.mean); u_mean = [u_stats.mean] 
        u_mean_err = (np.sqrt(1/np.sum(1/(lines['uerr'][~np.isnan(lines['uerr'])])**2)))
        uerr.append(u_mean_err); #u_mean_errlist = [u_mean_err]
        
        
        #p_mean = P_debi(q_mean, u_mean, q_mean_err, u_mean_err)
        #p_mean = np.sqrt(np.power(q_stats.mean, 2) + np.power(u_stats.mean, 2)) 
        #perr.append(P_err(p_mean, q_mean, u_mean, q_mean_err, u_mean_err))
        #perr.append((1/(np.power(p_mean, 2)))*(np.sqrt(np.power(q_stats.mean, 2)*np.power(q_mean_err, 2) + np.power(u_stats.mean, 2)*np.power(u_mean_err, 2))))
        
        P_OG = P_tot(q_mean, u_mean)
        Perr = P_err(P_OG, q_mean, u_mean, q_mean_err, u_mean_err)
        p_mean = P_debi3(P_OG, Perr)
        p.append(p_mean[0])
        perr.append(Perr[0])
        
    #return averaged in table form (for constistency with plotting functions)
    ave_data = Table([wave, flx, q, qerr, u, uerr, p, perr], names=['wave', 'flx', 'q', 'qerr', 'u', 'uerr', 'p', 'perr'])
    #ave_data = Table([wave, flx, q, qerr, u, uerr, p[0], perr[0]], names=['wave', 'flx', 'q', 'qerr', 'u', 'uerr', 'p', 'perr'])
    
    return(ave_data)

def get_weighted_mean(data_table):
    import numpy as np
    from astropy.table import Table
    from statsmodels.stats.weightstats import DescrStatsW
           
    wave = []; flx = []; q = []; qerr = []; u =[]; uerr =[]; p =[]; perr =[]; pa = []; paerr = [] #lists to store averages 
    #calculate averages **Note for list x, x[~numpy.isnan(x)] gets rid of nans in list and is used for each data list
    wave_ave = np.nanmean(data_table['wave'][~np.isnan(data_table['wave'])]); wave.append(wave_ave)
    flx_ave = np.nanmean(data_table['flx'][~np.isnan(data_table['flx'])]); flx.append(flx_ave)
    
    #for polarization calculate error-weighted mean and error on the wieghted mean 
    #Note: there is a list that holds all mean values and a list for q, u, qerr and uerr that updates with  
    #a single value each time b/c the function to calculate P and Perr take lists as inputs
    q_stats = DescrStatsW(data_table['q'][~np.isnan(data_table['q'])], weights=data_table['qerr'][~np.isnan(data_table['qerr'])], ddof=0)
    q.append(q_stats.mean); q_mean = [q_stats.mean] 
    q_mean_err = (np.sqrt(1/np.sum(1/(data_table['qerr'][~np.isnan(data_table['qerr'])])**2)))
    qerr.append(q_mean_err); #q_mean_errlist = [q_mean_err]
    
    u_stats = DescrStatsW(data_table['u'][~np.isnan(data_table['u'])], weights=data_table['uerr'][~np.isnan(data_table['uerr'])], ddof=0)
    u.append(u_stats.mean); u_mean = [u_stats.mean] 
    u_mean_err = (np.sqrt(1/np.sum(1/(data_table['uerr'][~np.isnan(data_table['uerr'])])**2)))
    uerr.append(u_mean_err); #u_mean_errlist = [u_mean_err]
    
    
    #p_mean = P_debi(q_mean, u_mean, q_mean_err, u_mean_err)
    #p_mean = np.sqrt(np.power(q_stats.mean, 2) + np.power(u_stats.mean, 2)) 
    #perr.append(P_err(p_mean, q_mean, u_mean, q_mean_err, u_mean_err))
    #perr.append((1/(np.power(p_mean, 2)))*(np.sqrt(np.power(q_stats.mean, 2)*np.power(q_mean_err, 2) + np.power(u_stats.mean, 2)*np.power(u_mean_err, 2))))
    
    P_OG = P_tot(q_mean, u_mean)
    Perr = P_err(P_OG, q_mean, u_mean, q_mean_err, u_mean_err)
    p_mean = P_debi3(P_OG, Perr)
    p.append(p_mean[0])
    perr.append(Perr[0])
    
    #return averaged in table form (for constistency with plotting functions)
    ave_data_table = Table([wave, flx, q, qerr, u, uerr, p, perr], names=['wave', 'flx', 'q', 'qerr', 'u', 'uerr', 'p', 'perr'])
    ave_data_table = get_pa_column(ave_data_table)
    
    return(ave_data_table)

#this is old and incorporated into get_max_min 6/21/24
def get_min(data_table, target_column, find_val = 0):
    """Function finds the minimum value in a column of a data table provided and returns corresponding wavelength [0]
    and index [1] to get whole row of data that pertains to the desired lowest value."""
    import numpy as np
    list_set = list(set(data_table[target_column]))
    list_set.sort()
    max_val = list_set[find_val]
    min_val = np.nanmin(data_table[target_column])
    for line in range(len(data_table)):
        if data_table[line][target_column] == min_val:
            min_wave = data_table[line]['wave']
            indx = line
            
            if target_column == 'p':
                pa = PA([data_table[line]['q']], [data_table[line]['u']])
                print("minimum Polarization: " + str(min_val)+ ", " + str(pa) + ", wavelength: " + str(min_wave))
    return(min_wave, indx, min_val)

#changed get_max to get_max_min 6/21/24 so if error thrown switch the function called 
def get_max_min(data_table, target_column, find_val = -1):
    """Function finds the maximum or minimum value in a column of a data table provided 
    returns corresponding wavelength [0], index [1] and value [2]. 
    Inputs: data_table - table of values, target_column - string of column name, find_val - float that specifies min or max 
        find_val = -1 specifies 1st largest value (by default), -2 specifies the second largest value and so on.
        find_val = 0 specifies 1st smallest value, 1 specifies second smallest value and so on. 
    To get whole row of data that pertains to the desired value use get_max_min_dat function."""
    import numpy as np
    length = len(data_table[target_column])
    list_dat = list(set(data_table[target_column]))
    list_dat.sort()
    
    if find_val < 0: #get max dat
        max_val = length + find_val
        val = list_dat[max_val]
    if find_val >= 0: #get min dat
        val = list_dat[find_val]
    #max_val = np.nanmax(data_table[target_column])
    for line in range(len(data_table)):
        if data_table[line][target_column] == val:
            val_wave = data_table[line]['wave']
            indx = line
            
            if target_column == 'p':
                pa = PA([data_table[line]['q']], [data_table[line]['u']])
                perr = P_err(val, [data_table[line]['q']], [data_table[line]['u']], [data_table[line]['qerr']], [data_table[line]['uerr']])
                #print("maximum Polarization: " + str(max_val)+ ", +/- " + str(perr) + ", wavelength: " + str(max_wave))
    return(val_wave, indx, val)

def get_max_min_dat(data, rest_wave, velocity_range, target_column = 'p', max_min = 'max', find_val = -1):
    from astropy.table import vstack, Table
    import numpy as np
    """ Returns data table for values correspond to the max/min of the specified target_column pertaining to the rest_wave and velocity range provided. 
        Each row in table are values for each table in data provided. Assumes data is provided as list of epoch data in order of days post max. 
        Default setting is getting the maximum value (max_min = 'max') of the total polarization (target_column = 'p')."""
    print(target_column + max_min +" values (corresponding to pmax) for rest line ", str(rest_wave))
    max_tab = Table()
    #max_tab = Table(names =  ['wave', 'vel', 'velerr', 'p', 'perr', 'pa', 'paerr', 'q', 'qerr', 'u', 'uerr', 'flx'])   

    dat_len = len(data)
    for ind, dat in enumerate(data):
        if dat_len > 1:
            print("**************epoch ", str(ind+1), "*****************************************")
        if velocity_range < 0: #blueside 
            min_wave = get_obswave(rest_wave, velocity_range)[0]
            line_reg = get_lines(dat, min_wave, rest_wave);# line_reg = get_pa_column(line_reg, adjust='on')
            line_reg = get_vel_column(line_reg, rest_wave)
        else: #redside
            max_wave, min_wave = get_obswave(rest_wave, velocity_range)
            line_reg = get_lines(dat, rest_wave, max_wave); #line_reg = get_pa_column(line_reg, adjust='on')
            line_reg = get_vel_column(line_reg, max_wave)
        #print(line_reg['wave'], line_reg['vel'])
        #print(line_reg['pa'], line_reg['paerr'])
        #print(line_reg['qerr'], line_reg['uerr'], line_reg['perr'])
        #print(line_reg['p'], line_reg['perr'])
        #print((line_reg['perr']/line_reg['p'])*(90/np.pi))
        
        binsize = line_reg['wave'][1] - line_reg['wave'][0]

        #max or min (first, second...) specified by find_val:
        vals = get_max_min(line_reg, target_column, find_val)
        data = dat.to_pandas() 
        i = vals[1] #idex of max/min value
        wave = np.round(vals[0], -1); max_tab['wave'] =[wave]
        vel = np.round(get_velocity(rest_wave, wave), -1); velerr = np.round(get_velocity(wave, wave + binsize/2), -1)
        max_tab['vel'] = [vel]; max_tab['velerr'] = [velerr]
        if 'flx' in data:
            print('HERE')
            flx = line_reg['flx'][i]
            max_tab['flx'] = [flx]
        if 'p' in data:
            p = np.round(vals[2], 2); perr = np.round(line_reg['perr'][i], 2)
            max_tab['p'] = [p]; max_tab['perr'] = [perr]
        if 'q' and 'qerr' and 'u' and 'uerr' in data:
            line_reg = get_pa_column(line_reg, adjust='on')
        if 'pa' in data:
            pa = np.round(line_reg['pa'][i], 1); paerr = np.round(line_reg['paerr'][i], 1)
            max_tab['pa'] = [pa]; max_tab['paerr'] = [paerr]
        if 'q' in data:
            q= np.round(line_reg['q'][i], 2); qerr = np.round(line_reg['qerr'][i], 2)
            max_tab['q'] = [q]; max_tab['qerr'] = [qerr]
        if 'u' in data:
            u = np.round(line_reg['u'][i], 2); uerr = np.round(line_reg['uerr'][i], 2)
            max_tab['u'] = [u]; max_tab['uerr'] = [uerr]
        #row = (wave, vel, velerr, p, perr, pa, paerr, q, qerr, u, uerr, flx)"""
        #max_tab.add_row(row)
        #print(max_tab[ind])    
    return(max_tab)

def get_norm_factor(norm_epoch, other_epoch, norm_region):
    """Calculates normalization factor for flux spectrum, returns this and the flux value it is normalized to. 
    Takes epoch data table to normalize to, epoch data table normalize, region of data to calculate max flx value for normalizing [min_wavelength, max_wavelength]
    ****NOTE: the two data tables input should be UNBINNED DATA****"""
    
    #get max flx value for epoch to normalize other epochs to 
    cont1 = get_lines(norm_epoch, norm_region[0], norm_region[1])
    max1 = get_max(cont1, 'flx') #get data line for max flux value 
    mf1 = cont1[max1[1]]['flx'] #get max flux value
    print("flux max to normalize to: " +str(mf1))

    cont2 = get_lines(other_epoch, norm_region[0], norm_region[1])
    max2 = get_max(cont2, 'flx')
    mf2 = cont2[max2[1]]['flx']
    print("flux max that needs normalizing: " + str(mf2))

    n12 = mf1/mf2 #calculate normalization factor 
    print("normalization factor: " + str(n12))
    print("check: " + str(mf2*n12))
    
    return(n12, mf1)

#get ellipse fit parameters for QU and plot ellipse 
def fit_ellipse(x, y):
    """

    Fit the coefficients a,b,c,d,e,f, representing an ellipse described by
    the formula F(x,y) = ax^2 + bxy + cy^2 + dx + ey + f = 0 to the provided
    arrays of data points x=[x1, x2, ..., xn] and y=[y1, y2, ..., yn].

    Based on the algorithm of Halir and Flusser, "Numerically stable direct
    least squares fitting of ellipses'.


    """

    D1 = np.vstack([x**2, x*y, y**2]).T
    D2 = np.vstack([x, y, np.ones(len(x))]).T
    S1 = D1.T @ D1
    S2 = D1.T @ D2
    S3 = D2.T @ D2
    T = -np.linalg.inv(S3) @ S2.T
    M = S1 + S2 @ T
    C = np.array(((0, 0, 2), (0, -1, 0), (2, 0, 0)), dtype=float)
    M = np.linalg.inv(C) @ M
    eigval, eigvec = np.linalg.eig(M)
    con = 4 * eigvec[0]* eigvec[2] - eigvec[1]**2
    ak = eigvec[:, np.nonzero(con > 0)[0]]
    return np.concatenate((ak, T @ ak)).ravel()


def cart_to_pol(coeffs):
    """

    Convert the cartesian conic coefficients, (a, b, c, d, e, f), to the
    ellipse parameters, where F(x, y) = ax^2 + bxy + cy^2 + dx + ey + f = 0.
    The returned parameters are x0, y0, ap, bp, e, phi, where (x0, y0) is the
    ellipse centre; (ap, bp) are the semi-major and semi-minor axes,
    respectively; e is the eccentricity; and phi is the rotation of the semi-
    major axis from the x-axis.

    """

    # We use the formulas from https://mathworld.wolfram.com/Ellipse.html
    # which assumes a cartesian form ax^2 + 2bxy + cy^2 + 2dx + 2fy + g = 0.
    # Therefore, rename and scale b, d and f appropriately.
    a = coeffs[0]
    b = coeffs[1] / 2
    c = coeffs[2]
    d = coeffs[3] / 2
    f = coeffs[4] / 2
    g = coeffs[5]

    den = b**2 - a*c
    if den > 0:
        raise ValueError('coeffs do not represent an ellipse: b^2 - 4ac must'
                         ' be negative!')

    # The location of the ellipse centre.
    x0, y0 = (c*d - b*f) / den, (a*f - b*d) / den

    num = 2 * (a*f**2 + c*d**2 + g*b**2 - 2*b*d*f - a*c*g)
    fac = np.sqrt((a - c)**2 + 4*b**2)
    # The semi-major and semi-minor axis lengths (these are not sorted).
    ap = np.sqrt(num / den / (fac - a - c))
    bp = np.sqrt(num / den / (-fac - a - c))

    # Sort the semi-major and semi-minor axis lengths but keep track of
    # the original relative magnitudes of width and height.
    width_gt_height = True
    if ap < bp:
        width_gt_height = False
        ap, bp = bp, ap

    # The eccentricity.
    r = (bp/ap)**2
    if r > 1:
        r = 1/r
    e = np.sqrt(1 - r)

    # The angle of anticlockwise rotation of the major-axis from x-axis.
    if b == 0:
        phi = 0 if a < c else np.pi/2
    else:
        phi = np.arctan((2.*b) / (a - c)) / 2
        if a > c:
            phi += np.pi/2
    if not width_gt_height:
        # Ensure that phi is the angle to rotate to the semi-major axis.
        phi += np.pi/2
    phi = phi % np.pi

    return x0, y0, ap, bp, e, phi


def get_ellipse_pts(params, npts=100, tmin=0, tmax=2*np.pi):
    """
    Return npts points on the ellipse described by the params = x0, y0, ap,
    bp, e, phi for values of the parametric variable t between tmin and tmax.

    """

    x0, y0, ap, bp, e, phi = params
    # A grid of the parametric variable, t.
    t = np.linspace(tmin, tmax, npts)
    x = x0 + ap * np.cos(t) * np.cos(phi) - bp * np.sin(t) * np.sin(phi)
    y = y0 + ap * np.cos(t) * np.sin(phi) + bp * np.sin(t) * np.cos(phi)
    return x, y

def mc_error_adjust(err_list, val_list, inst_err = 0.05):
    """ Uses a Monte Carlo method of generating a random number to determine 
    what value within a +/- error range to adjust an individual data point by.
    Takes in a list of error values and a corresponding list of point values ex. qerr, q. 
    Returns a list of error corrected point values. """
    from astropy.table import Table
    import random
    
    ecv = [] #error corrected values (q or u)
    for err, val in zip(err_list, val_list):
        #determine if intrinsic error is less than the instrumental error. 
        #Use instrumental error if so. 
        ###SNe meet 9/9/22 we decided that the instrumental error is not applicable for individual line regions.
        #if -inst_err < err < inst_err: 
            #err = inst_err
            
    #calculate randome error within error range using a random number     
        grn = random.uniform(0, 1) #generate random number
        srn = (grn*2) - 1 #shift random number to cover range -1 to +1 
        rev = srn*err #apply random number to error value to get random error value
        ecv.append(val + rev) #create list of error corrected values (OG point value + error)
        
    t = Table([ecv], names = ['ecv'])
    return(t['ecv'])

def get_loopy(Q, U, Qerr, Uerr, ax = 'None'):
    import matplotlib.pyplot as plt
    """Plots best fit ellipse given list of Q and U values and returns ellipse fit parameters
    (x0, y0) is the ellipse centre; 
    (ap, bp) are the semi-major and semi-minor axes, respectively; 
    e is the eccentricity; 
    phi is the rotation of the semi-major axis from the x-axis
    
    Uses functions get_ellipse, cart_to_pol and get_ellipse_pts from:
        https://scipython.com/blog/direct-linear-least-squares-fitting-of-an-ellipse/"""
    
    
    #shift all Q and U points to account for a radnomly determined error value within the points given error range 
    new_q = mc_error_adjust(Qerr, Q)
    new_u = mc_error_adjust(Uerr, U)
    
    if ax == 'None': #set axis if not given 
        ax = plt.gca()
    
    #using formula for ellipse: F(x,y) = ax^2 + bxy + cy^2 + dx + ey + f = 0
    coeffs = fit_ellipse(new_q, new_u) #get cartesian ellipse coefficiants (a, b, c, d, e, f) 
    #print('loop coordinates: a, b, c, d, e, f =', coeffs, 'from F(x,y) = ax^2 + bxy + cy^2 + dx + ey + f = 0')
    
    #convert cartesian ellipse coefficients to classic ellipse parameters:
    #where (x0, y0) is the ellipse centre; 
    #(ap, bp) are the semi-major and semi-minor axes, respectively; 
    #e is the eccentricity; 
    #phi is the rotation of the semi-major axis from the x-axis
    
    ellipse_parms = cart_to_pol(coeffs)
    print('x0, y0, ap, bp, e, phi =', ellipse_parms)

    x, y = get_ellipse_pts(cart_to_pol(coeffs)) #get x and y values to plot ellipse
    ax.plot(x, y, c = '#A65628', label = 'error corrected') #plot ellipse

    return(ellipse_parms)

def get_loopy_err(qs, us, qsig, usig):
    """ Runs get_loopy routine 100 times to calculate the average parameters output and the stdev of each. 
        Takes in lists of q, u, qsig and usig to fit the loop to. 
    """
    x0_parms = []; y0_parms = []; ap_parms = []; bp_parms = []; e_parms = []; phi_parms = []
    for i in range(100):
        parms = get_loopy(qs, us, qsig, usig)
        x0_parms.append(parms[0]); y0_parms.append(parms[1])
        ap_parms.append(parms[2]); bp_parms.append(parms[3])
        e_parms.append(parms[4]); phi_parms.append(parms[5])
    
    x0_ave = np.average(x0_parms); x0_sig = np.std(x0_parms)
    y0_ave = np.average(y0_parms); y0_sig = np.std(y0_parms)
    ap_ave = np.average(ap_parms); ap_sig = np.std(ap_parms)
    bp_ave = np.average(bp_parms); bp_sig = np.std(bp_parms)
    e_ave = np.average(e_parms); e_sig = np.std(e_parms)
    phi_ave = np.average(phi_parms); phi_sig = np.std(phi_parms)
    
    print(' ')
    print('Ran Loop-fitting routine 100x, returns [ave, stdev] for each parameter (in order):')
    print('[x0, y0, ap, bp, e, phi]')
    return([x0_ave, x0_sig], [y0_ave, y0_sig], [ap_ave, ap_sig], [bp_ave, bp_sig], [e_ave, e_sig], [phi_ave, phi_sig])

#find intersection of lower 5-sigma intersection 
def line_intersection(line1, line2):
    xdiff = (line1[0][0] - line1[1][0], line2[0][0] - line2[1][0])
    ydiff = (line1[0][1] - line1[1][1], line2[0][1] - line2[1][1])

    def det(a, b):
        return a[0] * b[1] - a[1] * b[0]

    div = det(xdiff, ydiff)
    if div == 0:
       raise Exception('lines do not intersect')

    d = (det(*line1), det(*line2))
    x = det(d, xdiff) / div
    y = det(d, ydiff) / div
    return x, y

def get_curvefit(x, y, xrange, sigma = 'none'):
    """Takes in X and Y data and returns xfit and yfit lists to plot an exponential curve fit
    as well as a and b parameters for the equation: y = ae^(bx) 
    The Third argument is an array sapnning the whole length in x you want to curve to be plotted (usually the same list as x)
    reference: https://rowannicholls.github.io/python/mathematics/curve_fitting/exponential.html"""
    if sigma != 'none':
        poly = np.polyfit(x, np.log(y), 1, w = 1/np.array(sigma), full = True, cov = 'unscaled') #polyfit weights by magnitude so to get smaller errors weighted more use standard 1/sigma_y**2
                                                                                    #cov = 'unscaled' returns covariance matrix unscaled b/c sigma are reliable estimates of the uncertainties 
    else: 
        poly = np.polyfit(x, np.log(y), 1, full = True, cov = True)
        #r = (np.polyval(np.polyfit(x, y, 1), x) - y) #residuals
        #gof = np.sum(r**2) #gof = rss = sum of the squares of the residuals (fit errors) https://stackoverflow.com/questions/15721053/whats-the-error-of-numpy-polyfit 
    
    a = np.exp(poly[0][1]); b = poly[0][0]
    x_fit = np.linspace(np.min(xrange), np.max(xrange), 100)
    y_fit = a*np.exp(b*x_fit)
    
    r = y - (a*np.exp(b*np.array(x)))
    chi2 = sum((r / sigma)**2) #unsure of results
    rss = np.sum(r**2) # sum of the squares of the residuals (fit errors) https://stackoverflow.com/questions/15721053/whats-the-error-of-numpy-polyfit 
    gof = rss
    
    return(x_fit, y_fit, a, b, gof)

def pow_bfit(x, y, xrange, y0 = 'none'):
    from scipy.optimize import curve_fit    
   
    def func_pow(x, a, b) :
        return a * np.power(x, b)
        
    if y0 == 'none':
        popt, pcov = curve_fit(func_pow, x, y) 
        x_fit = np.linspace(min(xrange), max(xrange)) #x pts
        fit = func_pow(x_fit, *popt) #y pts for fit
        r = y - func_pow(x, *popt) 
        chi2 = sum((r / y) ** 2) #not an accurate value
        rss = sum(r**2) #risidual sum of squares 
    else:
        popt, pcov = curve_fit(func_pow, x, y, sigma = y0, maxfev = 20000) #automatically optimizes fit from weights of errors https://stackoverflow.com/questions/27696324/using-scipy-optimize-curve-fit-with-weights
        x_fit = np.linspace(min(xrange), max(xrange)) #x pts
        r = y - func_pow(x, *popt)
        fit = func_pow(x_fit, *popt) #y pts for fit
        
        chi2 = sum((r / y0) ** 2)
        rss = sum(r**2)
    
    gof = rss #goodness of fit 
    return(x_fit, fit, gof, popt)


def lin_bfit(x, y, xrange = 'none', *args, **kwargs):
    #from scipy.odr import *
    w_x = kwargs.get('w_x')
    w_y = kwargs.get('w_y')
    #mydata = Data(x, y, wd=w_x, we=w_y) #we is the response variable, wd is the 
    mydata = RealData(x, y, sx=w_x, sy=w_y)
    #mydata = RealData(x, y, sx=1/np.power(w_x,2), sy=1/np.power(w_y,2))
    myodr = ODR(mydata, model=unilinear)
    #myodr.set_iprint(init=2, iter=2, final=2)
    myoutput = myodr.run()
    slope, intercept = myoutput.beta
    slope_err, intercept_err = myoutput.sd_beta
    if xrange != 'none':
        fit_x = np.linspace(xrange[0], xrange[1])
    else:
        fit_x = np.unique(x)
    fit_y = slope*fit_x + intercept
    #get chi**2 of fit
    chi_square = myoutput.res_var
    
    #print("fit parameter 1-sigma error")
    #print("Slope = "+ str(slope)+" +- "+str(slope_err))
    #print("Intercept = " + str(intercept) + " +- " + str(intercept_err))
    #print("———————————–")
    return(fit_x, fit_y, slope, intercept, slope_err, intercept_err, chi_square)

def get_slope_angle(slope, sigma_slope):
    import math
    """takes in slope and slope uncertainty as floats 
    retunrns angle of slope in degrees (divided by 2 for QU space)"""
    angle = (math.degrees(math.atan(slope)))/2
    sigma_angle = (math.degrees(sigma_slope/(1+(slope**2))))/2
    return(angle, sigma_angle)

def get_distance(x, y, x0, y0):
    """
    Return distance between point
    P[x0,y0] and a curve (x,y)
    """
    d_x = x - x0
    d_y = y - y0
    dis = np.sqrt( d_x**2 + d_y**2 )
    return(dis)

def get_min_distance(x, y, P, precision=5):
    """
    Compute minimum/a distance/s between
    a point P[x0,y0] and a curve given by (x,y)
    rounded at precision specified.
    
    ARGS:
        x, y      (array)
        P         (tuple)
        precision (int)
        
    Returns min indexes and distances array.
    To get min distance use the indexes to get that value from the distance array.
    """
    # compute distance
    d = get_distance(x, y, P[0], P[1])
    d = np.round(d, precision)
    # find the minima
    glob_min_idxs = np.argwhere(d==np.min(d)).ravel()
    return(glob_min_idxs, d)

#calculation of SN characteristic time
def solarMass_to_grams(solar_mass):
    return(solar_mass*(1.989*10**33))

def seconds_to_days(seconds):
    return(seconds/(60*60*24))

def sn_char_time(energy, solar_mass, particle_density = 1, k = 0.03):
    #Calculates the characteristic time of a SN in days
    #requires energy in ergs, ejecta mass in solar masses and particle number density 
    #(which is converted to proper g/cm**3 using the hydrogen atom mass)
    #in order to get the transition time between free-expansion to reverberation state (implosion time)
    #this answer should be multiplied by 2.4 
    #Reference Olmi B. 2023, Wheeler 2015 
    Hparticle_mass = 1.6735*10**-24 #grams
    mass = solarMass_to_grams(solar_mass) 
    density = particle_density*Hparticle_mass
    Con = 0.05
    #k = 0.03
    c_light = 3*10**10 #cm/s
    v = 1.25*10**9 #cm/s
    #char_t = seconds_to_days(energy**(-1/2) * mass**(5/6) * density**(-1/3)) #Olmi 2023 (Truelove 1999) method
    char_t = seconds_to_days((Con*k*mass**2/energy)**(1/2)) #Wheeler 2015 method
    return(char_t)

