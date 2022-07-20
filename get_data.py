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
    
    Perr = (1/(np.power(P, 2)))*(np.sqrt(np.power(Q, 2)*np.power(Qerr, 2) + np.power(U, 2)*np.power(Uerr, 2)))
    t = Table([Perr], names = ['perr'])
    return(t['perr'])

def P_debi(Q, U, Qerr, Uerr):
    from astropy.table import Table
    import numpy as np
    P = (np.sqrt(np.power(Q, 2) + np.power(U, 2) - 0.5*(np.power(Qerr, 2) + np.power(Uerr, 2))))
    t = Table([P], names = ['p'])
    return(t['p'])

def PA(Q, U, adjust_pa = 'off'):   
    from astropy.table import Table
    import numpy as np
    PA = np.degrees(.5*np.arctan2(U, Q)) #calculate PA
    
    if adjust_pa != 'on': #decide to add 180 to PA or not
        pass
    else:
        for i in range(len(PA)):
            if PA[i] < 0:
                PA[i] = PA[i]+180
            else:
                PA[i] = PA[i]+0
    t = Table([PA], names = ['PA'])
    return(t['PA'])

def PA_err(P, Q, U, Qerr, Uerr):
    """Calculates total Polarization Angle errors in quadrature from Q, Q errors and U, U errors values and returns errors as degrees.
    Data is returned as an astropy table column."""
    from astropy.table import Table
    import numpy as np
    PAerr = []
    for i in range(len(P)):
        err = 1/P[i]**2*(np.sqrt(Q[i]**2*Uerr[i]**2 + U[i]**2*Qerr[i]**2))
        PAerr.append(np.degrees(err))
    #PAerr = np.degrees((1/(np.power(P, 2)))*(np.sqrt(np.power(Q, 2)*np.power(Uerr, 2) + np.power(U, 2)*np.power(Qerr, 2))))
    
    t = Table([PAerr], names = ['pa_err'])
    return(t['pa_err'])

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
        output_data.write(dir+ new_filename+".txt", format = 'ascii')
    
    return(output_data)

def get_txtFITS(folder_path, flxfile, qfile, qerrfile, qsumfile, ufile, uerrfile, new_filename = "none"):    
    """
    **NOTE: this works for txt files made exactly from FITS files and uses the fits headings so if the file is not the txt version of the fits this will not work!**
    Takes in sepperate txt files for qsum, q, qerr, usum, u, uerr data (wavelength is the first column in all of them)
    and returns one table with all of it combined, if new_filename is provided this table is saved in a .txt file in the epoch folder
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
    wave = data1['wave']
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
 

def get_fits(file_path, names = ["wave", "q", "qerr", "qsum", "u", "uerr", "usum"], file_pattern = "\*.fits", new_filename = "none"):
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
    length = float(data[0].header['NAXIS1'])
    start = float(data[0].header['CRVAL1'])
    step = float(data[0].header['CDELT1'])
    stop = start + (length*step)
    waves = np.arange(start, stop, step);
    data_array = [waves] + data_array #add list of waves to begining of data array

    print(len(waves), len(data_array[1]))
    if len(waves) == len(data_array[1]):
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
    
    return((obs*conversion-rest*conversion)/(rest*conversion))*c

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
        qerr.append(q_mean_err); q_mean_errlist = [q_mean_err]
        
        u_stats = DescrStatsW(lines['u'][~np.isnan(lines['u'])], weights=lines['uerr'][~np.isnan(lines['uerr'])], ddof=0)
        u.append(u_stats.mean); u_mean = [u_stats.mean] 
        u_mean_err = (np.sqrt(1/np.sum(1/(lines['uerr'][~np.isnan(lines['uerr'])])**2)))
        uerr.append(u_mean_err); u_mean_errlist = [u_mean_err]
        
        #p_mean = P_tot(q_mean, u_mean)
        p_mean = np.sqrt(np.power(q_stats.mean, 2) + np.power(u_stats.mean, 2)) 
        p.append(p_mean)
        #perr.append(P_err(p_mean, q_mean, u_mean, q_mean_errlist, u_mean_errlist))
        perr.append((1/(np.power(p_mean, 2)))*(np.sqrt(np.power(q_stats.mean, 2)*np.power(q_mean_err, 2) + np.power(u_stats.mean, 2)*np.power(u_mean_err, 2))))
    
    #return averaged in table form (for constistency with plotting functions)
    ave_data = Table([wave, flx, q, qerr, u, uerr, p, perr], names=['wave', 'flx', 'q', 'qerr', 'u', 'uerr', 'p', 'perr'])
    
    return(ave_data)

def get_min(data_table, target_column):
    """Function finds the minimum value in a column of a data table provided and returns corresponding wavelength [0]
    and index [1] to get whole row of data that pertains to the desired lowest value."""
    import numpy as np
    min_val = np.min(data_table[target_column])
    for line in range(len(data_table)):
        if data_table[line][target_column] == min_val:
            min_wave = data_table[line]['wave']
            indx = line
            
            if target_column == 'p':
                pa = PA([data_table[line]['q']], [data_table[line]['u']])
                print("minimum Polarization: " + str(min_val)+ ", " + str(pa) + ", wavelength: " + str(min_wave))
    return(min_wave, indx)

def get_max(data_table, target_column):
    """Function finds the minimum value in a column of a data table provided and returns corresponding wavelength [0]
    and index [1] to get whole row of data that pertains to the desired lowest value."""
    import numpy as np
    max_val = np.max(data_table[target_column])
    for line in range(len(data_table)):
        if data_table[line][target_column] == max_val:
            max_wave = data_table[line]['wave']
            indx = line
            
            if target_column == 'p':
                pa = PA([data_table[line]['q']], [data_table[line]['u']])
                print("maximum Polarization: " + str(max_val)+ ", " + str(pa) + ", wavelength: " + str(max_wave))
    return(max_wave, indx)
