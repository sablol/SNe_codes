# -*- coding: utf-8 -*-
"""
Created on Fri Apr 23 17:42:06 2021

@author: sabri
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from pylab import *
from scipy import odr 
from scipy.optimize import curve_fit
from sklearn.linear_model import LinearRegression
import imageio
from get_data import * 
from scipy.stats import linregress
#from mpltools import color 

def plot_p_single(line_list, label_list, color_list, unbinned_t, bin_data_t, bin_size, flux_factor, title, save = "none", shift_list = "none"):
    """Creates a single panel plot of polarization and relative flux for one epoch. 
    Line, label and color lists are to specify features in spectra and may be given empty list []. 
    wl and flux are unbinned data list collected from get_data function.
    The binned lists can be gathered from get_binned. 
    Flux factor is a number to multiply the flux by so it is scaled to the polarization for presentation. 
    bin_size is a number provided for labeling. If a title is provided the figure will be saved as that. 
    If a shift_list (a list of wavelengths representing locations of shifted line features) is provided it will be plotted too.
    Title must be provided  and if save = on plot is saved under title in Single_epoch_tot_p folder.
    """
    
    #calculate total polarization
    p = bin_data_t['p']
    f = unbinned_t['flx']*flux_factor #scale flux 
    pmax = p.max()#max polarization value (for labels later)
    
    fig, ax = subplots(1, sharex = True, sharey = True, figsize = (14, 7))
    fig.suptitle(title, y = .95, fontsize = 22)
    #ax.axvspan(6462, 7102, color = 'lightgrey', label = '+/-15000km/s centered on 6800', alpha = .6) #alpha adjusts the transparency
    ax.step(bin_data_t['wave'], p, 'k', label = 'Polarization (Bin Size ' + str(bin_size) +"Å)")
    ax.fill_between(unbinned_t['wave'], f, color = 'powderblue', label = "Relative Flux") 
    ax.tick_params(axis='both', which='major', labelsize=14)
    ax.set_ylim([0, 5.2]); ax.set_xlim([4000, 7900])
    ax.legend(loc="upper left", fontsize = 14); 
    ax.set_xlabel("Wavelength(Å)", fontsize = 16); ax.set_ylabel("Percent Total Polarization", fontsize = 16)
    
    #ISP estimate areas
    #ax.axhline(y=.3, xmin=.335, xmax=.385, c= 'fuchsia')
    #ax.axhline(y=.1, xmin=.58, xmax=.62, c= 'fuchsia')
    
    for xc, c in zip(line_list, color_list):
        ax.axvline(x=xc, c=c) #plot veritcle line at location xc in linelist with color c in colorlist
    
    label_loc = 2 #initial label position    
    for xc, lab in zip(line_list, label_list):
        ax.text(xc+10, label_loc, lab, fontsize = 14, fontstyle = 'oblique', rotation = 'vertical') #add label at verticle position pmax next to line at xc, with label lab from label list 
        label_loc = label_loc + .2 #change label position each time so they don't overlap 
    
    #plot shifted lines in same color as rest lines if list is provided
    if shift_list != "none":
        for xc, c in zip(shift_list, color_list):
            ax.axvline(x=xc, c=c, linestyle = ':', linewidth = 3)
    if save != "none":
        fig.savefig("Single_epoch_tot_P/"+str(title))
    return(ax)


def plot_pannels(unbinned_t, bin_data_t, unbin_flux_factor, bin_flux_factor, adjust_pa = 'on', title = 'none'):
    
    """Creates the classic five plots (total polarization,polarized flux, %U, %Q, PA vs wavelength) with the flux 
    plotted in a shaded region of each pannel. 
    
    Takes in two astropy tables for unbinned and binned data
    unbinned_file must contain a column call 'flx' (so make sure to pass that to get_fits if that's used)
    
    Flux factors are what to multiply the unbinned and binned flux by to be on the same scale as the polarization. 
    adjust_pa is automatically set to add 180 for values less than 0, if you don't want this adjust_pa = 'off'.
    
    Provide a title to save plot otherwise the plot will be titled 'none' and will not be saved. 
    """
    
    import numpy as np
    import matplotlib.pyplot as plt
    
    Q = np.array(bin_data_t['q'])*100 #convert fraction to percent 
    U = np.array(bin_data_t['u'])*100
    
    F = np.array(unbinned_t['flx'])*unbin_flux_factor #Adjust scale of relative flux to magnitude of polarizations
    
    P_Flux = bin_data_t['p']*bin_data_t['flx']*bin_flux_factor
    
    PA = np.degrees(.5*np.arctan2(U, Q)) #calculate PA
    
    if adjust_pa != 'on': #decide to add 180 to PA or not
        pass
    else:
        for i in range(len(PA)):
            if PA[i] < 0:
                PA[i] = PA[i]+180
            else:
                PA[i] = PA[i]+0

    PA_ave = [] #creat list for plotting line of average PA 
    for i in range(len(PA)): PA_ave.append(np.nanmean(PA)) #calculate average PA
    
    fig, (ax0, ax1, ax2, ax3, ax4) = plt.subplots(5, sharex = True, figsize=(20,14))
    fig.suptitle(title, fontsize=24, y=.95)
    ax0.step(bin_data_t['wave'], bin_data_t['p'], 'r'); ax0.fill_between(unbinned_t['wave'], F, color = 'skyblue', label = 'Relative Flux'); ax0.set_ylabel("%Polarization", fontsize = 'medium'); ax0.legend(fontsize='medium')
    ax1.step(bin_data_t['wave'], P_Flux, 'r'); ax1.fill_between(unbinned_t['wave'], F, color = 'skyblue'); ax1.set_ylabel("Polarized Flux", fontsize = 'medium')
    ax2.step(bin_data_t['wave'], Q, 'r'); ax2.fill_between(unbinned_t['wave'], F, color = 'skyblue'); ax2.set_ylabel("% Q", fontsize = 'medium')   
    ax3.step(bin_data_t['wave'], U, 'r'); ax3.fill_between(unbinned_t['wave'], F, color = 'skyblue' ); ax3.set_ylabel("% U", fontsize = 'medium')   
    ax4.step(bin_data_t['wave'], PA, 'r'); ax4.plot(bin_data_t['wave'], PA_ave, 'k', label = 'Average PA: '+ str(round(np.nanmean(PA))) + u'\xb0', linewidth = 3);ax4.set_ylabel("Position Angle", fontsize = 'medium'); ax4.legend(fontsize = "medium"); ax4.set_xlabel("Wavelength(Å)")
    
    plt.subplots_adjust(wspace=0, hspace=0)
    #if title != 'none':
        #fig.savefig("Pannel_Plots/"+ title)
        #print("Plot has been saved under filename " + str(title) + " in the Pannel_Plots directory")
       
    return (print("provide title to save plot"))

def best_fit(Q, U, Qsig, Usig):
    from scipy import odr 
    from scipy.optimize import curve_fit
    from sklearn.linear_model import LinearRegression
    
    """Creates an error-weighted best fit line for Q/U plots given 4 lists of percent. 
    """
    
    def func(b, x):   
        # Linear function y = m*x + b
        # b is a vector of the parameters.
        # x is an array of the current x values.
        return b[0] * x + b[1]

    #model object
    lin_model = odr.Model(func)

    #creat real data object
    data = odr.RealData(Q, U, sx=Qsig, sy=Usig)

    # Set up ODR with the model and data.
    odr = odr.ODR(data, lin_model, beta0=[0., 1.])

    # Run the regression.
    out = odr.run()

    #fit peramaters
    popt = out.beta
    perr = out.sd_beta
    print("fit parameter 1-sigma error")
    print("Slope = "+ str(popt[0])+" +- "+str(perr[0]))
    print("Intercept = " + str(popt[1]) + " +- " + str(perr[1]))
    print("———————————–")
    x_fit = np.linspace(min(Q), max(Q))
    fit = func(popt, x_fit)
    
    #links for weighted fit:
    #https://micropore.wordpress.com/2017/02/07/python-fit-with-error-on-both-axis/
    #https://docs.scipy.org/doc/scipy/reference/odr.html
    
    return (x_fit, fit)

def plot_QU(binned_data, title, size = 3.2, save = 'none', epoch_labels = 'none'):
    
    """Creates a Q vs U plot with points color coded by wavelength and a line of best fit 
    that includes errors when epoch_labels = 'none' (most common use and is set by defualt). 
    
    Takes in an astropy Table with 'wave', 'q', 'u', 'qerr', 'uerr' columns. 
    size is auto set to .05 otherwise change for tighter fit. 
    Title must be provided  and if save = on the plot is saved under that in a QU_Plots directory.  
    A best fit including weighed errors is calculated from lists provided using the function best_fit.
    
    Optionally one can plot Q/U values without specifying wavelength by color if epoch_labels is set equal to something. 
    For example to plot the single average Q/U point for three epochs on the same plot 
    epoch_labels = ['epoch_1', 'epoch_2', 'epoch_3'] and each point will be labeled accordingly. 
    """
    
    wl = binned_data['wave']
    #convert fraction to percentage and ignore nan values
    Q = [x for x in binned_data['q'] if np.isnan(x) == False] 
    U = [x for x in binned_data['u'] if np.isnan(x) == False]
    Qsig = [x for x in binned_data['qerr'] if np.isnan(x) == False]
    Usig = [x for x in binned_data['uerr'] if np.isnan(x) == False]
    
    wl, Q = zip(*zip(wl, Q)) #makes sure wavelength and other lists are the same length for plotting
                             #by pairing values up in each then dropping unpaired values 
    
    bfit = best_fit(Q, U, Qsig, Usig)
    x_fit = bfit[0]
    fit = bfit[1]
    
    #calculate average error
    Qave_err = np.nanmean(Qsig)
    Uave_err = np.nanmean(Usig)

    #plot
    fig, ax = plt.subplots(1, figsize=(10,10))
    rcParams["font.size"]= 20
    fig.suptitle(title, fontsize=20, y=.85, x=.44) 
    ax.axhline(0, color = 'k', linewidth = 1, zorder = 1)
    ax.axvline(0, color = 'k', linewidth=1, zorder = 2)
    ax.plot(Q, U, c = 'lightgrey', linewidth = 1, zorder = 3 )
    ax.axis('square'); ax.set_xlim([-float(size), float(size)]); ax.set_ylim([-float(size), float(size)]) #comment out to see zoomed in 
    ax.set_xlabel("% Q Polarization", fontsize=18); ax.set_ylabel("% U Polarization", fontsize=18)

     
    if epoch_labels == 'none':
        #if plotting a a single epoch of data plot colors accoring to wavelength
        im = ax.scatter(Q,U, marker = 'o', c=wl, s=20, zorder = 4) #c=color set to a different color for each point in wavelength array, s=size 
        fig.colorbar(im).ax.set_ylabel('Wavelength(Å)', fontsize=20) #shows colorbar and labels it
        ax.plot(x_fit, fit, "r", lw=2, label="Best Fit", zorder = 5) #best fit line
        ax.errorbar((-float(size) + .005), (float(size) - .0058), xerr = Qave_err, yerr = Uave_err) #add error bar example in location based on grid size 
        #ax.text((-float(size)+ .008), (float(size) - .0055), "Average Error", fontsize = 13)
    else:
        #if ploting average data points (wavelength doesn't matter) plot all points as same color 
        im = ax.scatter(Q,U, marker = 'o', s=20, zorder = 5)
        ax.errorbar(Q, U, yerr=Usig, xerr=Qsig, ecolor="grey", hold=True, fmt="none", zorder = 4 ) #comment out to see points by wavelength
        for i, txt in enumerate(epoch_labels):
            ax.annotate(txt, xy = (Q[i]+0.0005, U[i]), fontsize = 14) #label each data point according to epoch_labels provided     
    
    ax.legend(loc="lower left",fontsize='x-small')
    if save != 'none':
        fig.savefig("QU_Plots/"+ title)
        print("figure saved in QU_Plots directory under" + str(title))
    else:
        print("set save = on to save figure")
    fig.show()
     
    #return (print("Q vs U plot has been saved under filename " + str(title) + " in the QU_Plots directory"))
    return() 
   
def plot_line_velspace(data_t, rest_wave, velocity, flx_factor, epoch, save = "off", plot = "off"):
    """Creates a velocity space plot centered on desired line. 
    Takes in a data table (from Bin_data or get_txtFITS), the wavelength of the line to center plot on, 
    the velocity range to go out to on either end, the flux factor so flux fits to scale with polarization, 
    an epoch for the plot title (no spaces) and option to save if save = on to folder VelSpace_plots/epoch+rest_wave
    (save is automatically off so plot will not save)"""

    title = str(epoch) + " Velocity Region Centered on " + str(rest_wave)
    vel=[] #creat list for velocity values
    wave_range = get_obswave(rest_wave, velocity) #calculate wavelengths on either side of rest that pertain to velocity region desired
    line_data = get_lines(data_t, wave_range[1], wave_range[0]) #get lines from data table in between wavelengths in velocity range
    #print(line_data)
    for i in line_data['wave']:
        vel.append(get_velocity(rest_wave, i)) #calculate velocities for each wavelength in selected data region
    
    if plot != "off":
        #plot total polarization and flux in velocity range
        plt.figure(figsize= (10,7)); plt.title(title, fontsize = '20')
        plt.plot(vel, line_data['p'], color = 'k', label = "Total Polarization")
        plt.fill_between(vel, line_data['flx']*flx_factor, color = 'skyblue', label = 'Relative Flux')
        plt.ylabel('Percent Total Polarization', fontsize = '20'); plt.xlabel('Velocity (km/s)',fontsize = '20'); plt.legend(fontsize = '14')
    
    if save != "off":
        #plt.savefig("VelSpace_plots/"+str(epoch)+str(rest_wave)) #option to save plot
        plt.savefig("VelSpace_plots/"+str(epoch)+str(rest_wave))
    
    return(line_data)

def plot_all_data(unbin_epoch_tables_list, epoch_tables_list, flux_adjust_list, pflux_adjust_list, epoch_labels_list, title = "none", xmin = "none", xmax = "none",):
    #find xmin, xmax from first epoch values if not otherwise provided 
    if xmin == 'none':
        xmin = np.min(epoch_tables_list[0]['wave'])
    if xmax == 'none':
        xmax = np.nanmax(epoch_tables_list[0]['wave'])
        
    #plot polarization and flux for each epoch stacked 
    fig, (ax0, ax1, ax2, ax3, ax4, ax5) = plt.subplots(6, sharex = True, figsize = (20, 24))
    
    if title != 'none':
        fig.set_title(title, fontsize = 16)
    #fig.set_xlim([xmin, xmax]) 
    ax0.set_ylabel('Relative Flux'); ax1.set_ylabel('% Polarization'); ax2.set_ylabel('polarized flux'); ax3.set_ylabel('% Q'); ax4.set_ylabel('% U'); ax5.set_ylabel('PA')
    ax5.set_xlabel('Wavelength  (Å)')  

    count = 0
    #plot flux with proper adjustment to be on same scale
    for epoch in unbin_epoch_tables_list:  
         color = iter(cm.rainbow(np.linspace(0, 1, count))); 
         ax0.plot(epoch['wave'], epoch['flx']*flux_adjust_list[count], label = epoch_labels_list[count], lw = 1)
         ax0.legend()
         count = count + 1
        
    #plot polarization data 
    count = 0
    for epoch in epoch_tables_list:
        ax1.step(epoch['wave'], epoch['p'])
        ax2.step(epoch['wave'], epoch['p']*epoch['flx']*pflux_adjust_list[count]); count = count+1
        ax3.step(epoch['wave'], epoch['q']) #already in %p 
        ax4.step(epoch['wave'], epoch['u'])
        ax5.step(epoch['wave'], PA(epoch['q'], epoch['u'])) #PA function found in get_data.py
     
    return(fig)

def plot_all_flx(unbin_epoch_tables_list, flx_adjust_list, epoch_names, num_epochs):
    
    plt.figure(figsize=(15, 12))
    plt.subplots_adjust(hspace=0)
    plt.rc('axes', labelsize= 12)    # fontsize of the x and y labels
    plt.rc('xtick', labelsize= 12)    # fontsize of the tick labels
    plt.rc('ytick', labelsize= 12)    # fontsize of the tick labels
    plt.rc('legend', fontsize= 12)    # legend fontsize
    #plt.xlabel('Wavelength  (Å)'), plt.ylabel('Relative Flux')
    # loop through the length of list of unbinned epoch data tables and keep track of index
    for n, epoch in enumerate(unbin_epoch_tables_list):
        # add a new subplot iteratively
        ax = plt.subplot(num_epochs, 1, n+1)#stack them
        ax.set_xlim(4000, 8000); ax.set_xticks([])#share same x-axis scale but don't show for each of them 
        ax.set_xlabel('Wavelength  (Å)', fontsize = 12)
        ax.plot(epoch['wave'], epoch['flx']*flx_adjust_list[n], c = 'k', label = str(epoch_names[n]) ) #label each one by epoch
        ax.set_yticks([1,3,5])
        #ax.axvline(4924, c = 'b')
        plt.legend() 
        if n == 0: 
            ax.annotate('He I', (5600, 5.85))
            ax.annotate('Fe II', (4750, 4.25))
            ax.annotate('$\\bigoplus$', (6850, 4.25))
            ax.annotate('$\\bigoplus$', (7600, 2))
        if n == 5: 
            ax.annotate('Mg I', (4571, 3.25))
            ax.annotate('Fe II (λλλ)', (4900, 3.25))
            ax.annotate('O I', (5577, 2.5))
            ax.annotate('NaID', (5890, 2.5))
            ax.annotate('O I', (6364, 4.5))
        if n == len(epoch_names)/2-1:   
           ax.set_ylabel('Relative Flux', fontsize = 12) 
    ax.set_xticks(np.arange(4000, 8000, 250))#show only x-axis for last plot
    return()

def plot_all_epochs(epoch_tables_list, flx_adjust_list, epoch_labels_list, title = "none", xmin = "none", xmax = "none", regions_list = 'none', line_list = 'none', label_list = 'none', color_list='none', polflux = 'off'):
    """Creates a plot of all epochs total polarization and flux stacked (spaced 5% appart). 
    Requires list of epoch data tables (outputs of Bin_data), 
    list of flux adjustment values so they are scaled to polarization levels, 
    list of epoch labels as strings. 
    Optional arguments are title, xmin (that epoch labels are located at), xmax and 
    list of regions of interest (like continuum) where each region is a list of two numbers (the bounds of the region).
    ex: regions_list = [[5000, 5500], [7000, 7300]]
    There is currently no save option so plots must be save from output.
    The function returns the axis so anything can be added to the plot after by:
        image = plot_all_epochs(blah, blah....)
        image.axvline(x=...)
    """
    
    #count number of epochs
    num_epochs = 0
    for count, item in enumerate(epoch_tables_list):
        num_epochs = num_epochs +1
    
    #find xmin, xmax from first epoch values if not otherwise provided 
    if xmin == 'none':
        xmin = np.min(epoch_tables_list[0]['wave'])
    if xmax == 'none':
        xmax = np.nanmax(epoch_tables_list[0]['wave'])
    
    #plot polarization and flux for each epoch stacked 
    fig, (ax) = plt.subplots(1, sharex = True, figsize = (10, 10))
    
    if title != 'none':
        ax.set_title(title, fontsize = 12)
    #ax.set_xlim([xmin, xmax])                                          
    count = 0 #count each epoch to index correct flx adjust and epoch label 
    fmax_list = [] #creat list to store maximum values to position labels

    if polflux != 'off':
        space = (num_epochs-1)*4 #calculate total space needed to plot epochs stacked
        ax.set_xlabel("Wavelength  (Å)", fontsize = 14); ax.set_ylabel("Relative Polarized Flux", fontsize = 14)
        for epoch in epoch_tables_list:
            P_Flux = epoch['p']*epoch['flx']*flx_adjust_list[count]
            ax.step(epoch['wave'], P_Flux + space)
            #ax.annotate(epoch_labels_list[count], xy = (xmin + 50, space +2), fontsize = 14) 
            count = count +1; space = space - 5

    else:
        space = (num_epochs-1)*5 #calculate total space needed to plot epochs 5% appart
        ax.set_xlabel("Wavelength  (Å)", fontsize = 12); ax.set_ylabel("Percent Polarization", fontsize = 12)
        for epoch in epoch_tables_list:
            p = ax.step(epoch['wave'], epoch['p'] + space, label = epoch_labels_list[count] + " (P+"+str(np.round(space))+"%)") #plot polarization of epochs spaced 5% appart       
            ax.plot(epoch['wave'], epoch['flx']*flx_adjust_list[count] + space, color = 'k', linestyle = ':')#plot flux with adjusted scale and inline with polarization
            fmax = np.nanmax(epoch['flx'])*flx_adjust_list[count] + space; fmax_list.append(fmax) #find max value and store
            #ax.legend(loc = 'upper right')
            ax.annotate(epoch_labels_list[count] + " \n (P+"+str(np.round(space))+"%)", xy = (xmin + 50, space+3 ), fontsize = 10, c=p[0].get_c()) #label epoch and polarization factor (color code label to line color)
            count = count +1; space = space - 5

        #highlight special regions (if lists of region bounds are provided)
            if regions_list != 'none':
                count_lists = 0 #index first list in list of lists

                for region in regions_list:
                    aves = get_aves([epoch], region) #calculate averages for data in region 
                    ave_p = str(np.round(aves['p'][0], 3)); #get average total polarization
                    print(ave_p)
                    ave_perr = str(np.round(aves['perr'][0], 4)) #get average total polarization error 
                    ax.axvspan(region[0], region[1], color = 'lightgrey') #highlight region 
                    #ax.annotate(ave_p + "% \n \u00B1" + ave_perr, xy =(region[0] + 50, space+8 ), fontsize = 9) #lable polarization and error in region where xy = is a really complicated way of locating the label        
                count_lists = count_lists + 1 #move to next region list index
            
    #add lines for features of interest 
    if line_list != 'none':
        for xc, c in zip(line_list, color_list):
            ax.axvline(x=xc, c=c) #plot veritcle line at location xc in linelist with color c in colorlist
    
        label_loc = 3.5 #initial label position    
        for xc, lab in zip(line_list, label_list):
            ax.text(xc-30, 28, lab, fontsize = 12, rotation = 'vertical') #add label at verticle position pmax next to line at xc, with label lab from label list 
            #ax.text(xc+10, fmax_list[0] - .1, lab, fontsize = 14) #rotation = 'vertical') #add label at verticle position pmax next to line at xc, with label lab from label list 
            #label_loc = label_loc + .2 #change label position each time so they don't overlap 

        count = count + 1 #move to next epoch index 
        space = space - 5 #reduce space so that oldest epic is on bottom of plot
        
    return(ax) 

def QU(binned_data, ax = None, bfit = 'none', size = 4.5, color_data = 'wl', epoch_labels = 'none', cmap_choice = 'turbo', bfit_c = 'grey'):
    """Returns QU plot given a binned_data table. The return object should be added to a figure in order
    to do fig.colorbar(QU_return) and fig.suptitle to make it look pretty. Multiple QU plots can be added 
    as subplots of a figure by: 
        import matplotlib.pyplot as plt
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize = (12, 5))
        e1 = QU(epoch1_20, 'E1', ax1)
        e2 = QU(epoch2_20, 'E2', ax2)
        fig.colorbar(e1).set_label('Wavelength(Å)') 
    A best fit line can be calculated outside the function and passed to it to exclude points that are in the QU plot. 
    If bfit is not specified the bestfit in calculated in the function using all the QU points. """
    
    if ax == None:
        ax = plt.gca()
    
    wl = binned_data['wave']
    #convert fraction to percentage and ignore nan values
    Q = [x for x in binned_data['q'] if np.isnan(x) == False] 
    U = [x for x in binned_data['u'] if np.isnan(x) == False]
    Qsig = [x for x in binned_data['qerr'] if np.isnan(x) == False]
    Usig = [x for x in binned_data['uerr'] if np.isnan(x) == False]
    
    wl, Q = zip(*zip(wl, Q)) #makes sure wavelength and other lists are the same length for plotting
                             #by pairing values up in each then dropping unpaired values 
    
    if bfit == 'none':
        bfit = best_fit(Q, U, Qsig, Usig)
        
    x_fit = bfit[0]
    fit = bfit[1]
    
    #calculate average error
    Qave_err = np.nanmean(Qsig)
    Uave_err = np.nanmean(Usig)
    print(Qave_err, Uave_err)
    if Qave_err < 0.1:
        Qave_err = 0.1
    if Uave_err < 0.1:
        Uave_err = 0.1 

    #plot
    ax.axhline(0, color = 'k', linewidth = 1, zorder = 1)
    ax.axvline(0, color = 'k', linewidth=1, zorder = 2)
    ax.plot(Q, U, c = 'lightgrey', linewidth = 2, zorder = 3 )
    ax.axis('square'); ax.set_xlim([-float(size), float(size)]); ax.set_ylim([-float(size), float(size)]) #comment out to see zoomed in 
    #ax.set_xlabel("q(%)", fontsize = 18); ax.set_ylabel("u(%)", fontsize = 18); ax.yaxis.set_label_coords(-.13, .5)
    #ax.tick_params(labelsize=18)
    
    if epoch_labels == 'none':
        #if plotting a a single epoch of data plot colors accoring to wavelength
        im = ax.scatter(Q,U, marker = 'o', c=wl, s=20, zorder = 4, cmap = cmap_choice) #c=color set to a different color for each point in wavelength array, s=size 
        #if another column of data is provided for colorcoding use that (for example 'vel' column from get_vel_column to colorcode in velocity space)
        if color_data != 'wl':
            im = ax.scatter(Q,U, marker = 'o', c=color_data, s=20, zorder = 4, cmap = cmap_choice) #c=color set to a different color for each point in wavelength array, s=size 
        #fig.colorbar(im).ax.set_ylabel('Wavelength(Å)', fontsize=12) #shows colorbar and labels it
        #ax.plot(x_fit, fit, bfit_c, lw=2, label="Best Fit", zorder = 5) #best fit line
        ax.errorbar((-float(size) + .5), (float(size) - .58), xerr = Qave_err, yerr = Uave_err) #add error bar example in location based on grid size 
        ax.text((-float(size)+ .8), (float(size) - .55), "Average Error", fontsize = 10)
        #ax.legend(loc="lower left", fontsize='medium')
    else:
        #if ploting average data points (wavelength doesn't matter) plot all points as same color 
        im = ax.scatter(Q,U, marker = 'o', s=20, zorder = 5)
        ax.errorbar(Q, U, yerr=Usig, xerr=Qsig, ecolor="grey", hold=True, fmt="none", zorder = 4 ) #comment out to see points by wavelength
        for i, txt in enumerate(epoch_labels):
            ax.annotate(txt, xy = (Q[i]+0.0005, U[i]), fontsize = 14) #label each data point according to epoch_labels provided     
    
    return(im)

def plot_QU_spec_comp2(unbindata, bindata, flx_adjust, region_lists, region_names, comp_color = 'whole_region'):
    """Creates a plot with the flux and polarization spectrum of a single epoch of data as the bottom pannel
    and two square plots comparing Q/U data of selected 2 regions in that same epoch above (hence comp2). 
    Q/U plots are color coded in velocity space with the colorbar dictated by the second region given. 
    Inputs required are: unbinned data table from get_fits, bindata table from Bin_data, 
    flx_adjust = a single number for flx adjust to adjust flx spectrum to polarization level for bottom pannel, 
    region_lists = a list of two lists for the two line regions to compare where in each list the entries are 
    [blue_limit_wave, rest_wave, red_limit_wave], region_names = a list of region names as strings and the region, 
    comp_color = a string value (blue or red) that specifies what side of the rest wavelength you want to compare in Q/U space.
    If no color is given the default is to compare the whole line region from the blue limit to the red limit
    
    ***NOTE: whatever region has complete data shoud be the second region listed, 
    if the second region does not have data spanning the whole region specified the Q/U velocity color scale will be inaccurate.  
    """
    
    fig = plt.figure()
    fig.set_figheight(8)
    fig.set_figwidth(10)

    ax1 = plt.subplot2grid(shape=(2,2), loc=(1,0), colspan = 2)
    ax2 = plt.subplot2grid(shape=(2,2), loc=(0,0), colspan = 1)
    ax3 = plt.subplot2grid(shape=(2,2), loc=(0,1), colspan = 1)
    ax_list = [ax2, ax3] #put in list for QU plots loop
    
    #bottom panel - spectra 
    ax1.fill_between(unbindata['wave'], unbindata['flx']*flx_adjust, color = 'silver')
    ax1.step(bindata['wave'], bindata['p'], c = 'k')
    ax1.set_ylim([0, 5.2])
    ax1.set_ylabel('Percent Total Polarization')
    ax1.set_xlabel('Wavelength (Å)')
    
    #top sqaures - line region QUs
    ax2.annotate(region_names[0], (3, 4))
    ax2.set_ylabel("%U")
    ax2.set_xlabel("%Q", loc = 'right')
    ax2.add_patch(plt.Circle((0, 0), radius = .32, edgecolor='red', facecolor = 'none'))#isp

    ax3.annotate(region_names[1], (3, 4))
    ax3.add_patch(plt.Circle((0, 0), radius = .32, edgecolor='red', facecolor = 'none'))#isp

     
    for region, ax in zip(region_lists, ax_list):
        ax1.axvspan(region[0], region[1], color = 'paleturquoise', alpha = .5)
        ax1.axvspan(region[1], region[2], color = 'lightcoral', alpha = .4)
        
        #get QU data in velocity space for blue region of lines
        if comp_color == 'blue':
            reg = get_lines(bindata, region[0], region[1])
            reg_vel = get_vel_column(reg, region[1])
            reg_qu = QU(reg_vel, ax = ax, color_data = reg_vel['vel'], cmap_choice = 'winter')
            
        #get QU data in velocity space for red region of lines
        if comp_color == 'red':
            reg = get_lines(bindata, region[1], region[2])
            reg_vel = get_vel_column(reg, region[1])
            reg_qu = QU(reg_vel, ax = ax, color_data = reg_vel['vel'], cmap_choice = 'spring')
            
        if comp_color == 'whole_region':
            reg = get_lines(bindata, region[0], region[2])
            reg_vel = get_vel_column(reg, region[1])
            reg_qu = QU(reg_vel, ax = ax, color_data = reg_vel['vel'], cmap_choice = 'turbo')
            
    plt.subplots_adjust(left=0, bottom=0, right=.85, top=.99, wspace=0.01, hspace=0.13)
    cbar_ax = fig.add_axes([.85, 0.5, 0.02, 0.5])
    fig.colorbar(reg_qu, cax=cbar_ax).set_label('Velocity(km/s)', rotation= 90) 
    return()

def make_gif(file_list, movie_name, slide_time):
    """makes a gif out of files provided in the file_list. 
    the Gif is saved under the movie_name.GIF (make sure to include .GIF at the end of the movie name) 
    and the slide_time specifies how long to spend on each image"""
    
    images = []

    for filename in file_list:
        images.append(imageio.imread(filename))
    imageio.mimsave(str(movie_name), images, duration = slide_time)
    return(print("GIF saved"))

