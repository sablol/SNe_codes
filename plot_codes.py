# -*- coding: utf-8 -*-
"""
Created on Fri Apr 23 17:42:06 2021

@author: sabri
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.collections import LineCollection
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)
#from matplotlib.colors import DivergingNorm
import matplotlib.colors as mcolors
from pylab import *
from scipy import odr 
from scipy.optimize import curve_fit
from sklearn.linear_model import LinearRegression
import imageio
from get_data import * 
from scipy.stats import linregress
from astropy.table import vstack, Table
#from mpltools import color 

def plot_p_single(line_list, label_list, color_list, unbinned_t, bin_data_t, bin_size, flux_factor, title, plot = 'p', save = "none", shift_list = "none"):
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
    ax.tick_params(axis='both', which='major', labelsize=14)
    ax.set_xlim([4000, 8000])
    ax.legend(loc="upper left", fontsize = 14); 
    ax.set_xlabel("Wavelength(Å)", fontsize = 16)
    ax.plot(unbinned_t['wave'], f, c = 'k')
    ax.fill_between(unbinned_t['wave'], f, color = 'powderblue', label = "Relative Flux") 
    if plot == 'pa':
        ax.set_ylim([-100, 100]);ax.set_ylabel("Polarization Angle", fontsize = 16)
        PAs = PA(bin_data_t['q'], bin_data_t['u'])
        PAerrs = PA_err(p, bin_data_t['q'], bin_data_t['u'], bin_data_t['qerr'], bin_data_t['uerr'])
        ax.step(bin_data_t['wave'], PAs, 'k', label = 'Polarization Angle (Bin Size ' + str(bin_size) +"Å)")   
    else:
        ax.set_ylim([0, 4]); ax.set_ylabel("Percent Total Polarization", fontsize = 16)
        ax.step(bin_data_t['wave'], p, 'k', label = 'Polarization (Bin Size ' + str(bin_size) +"Å)")
        
    #ISP estimate areas
    #ax.axhline(y=.3, xmin=.335, xmax=.385, c= 'fuchsia')
    #ax.axhline(y=.1, xmin=.58, xmax=.62, c= 'fuchsia')
    
    for xc, c in zip(line_list, color_list):
        ax.axvline(x=xc, c=c) #plot veritcle line at location xc in linelist with color c in colorlist
    
    label_loc = 3 #initial label position    
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

def plot_spec_shift(unbin_data, data, velocity, binsize, flx_adjust, title, shift_spec = 'flux'):
    """Plots chosen (%p by defualt or flux if specified) spectrum shifted by given velocity (- for blueshift, + for redshift).
    Takes unbinned data table (and recalculated wavlengths) to plot shifted flux spectrum. 
    Binned data for %p to be plotted at rest.
    Velocity to shift by. 
    Bin size of data. 
    Amount to multiply flux by to fit it on plot.
    Title of plot.  
    #Option to plot shifted %p spectrum instead of flux spectrum if shift_spec == 'p'
    #Option to plot %polflux calculated from rest the shifted instead of flux spectrum if shift_spec == 'polflux'
    
    Calls plot_p_single to creat plot and retruns this figure."""
    if shift_spec == 'flux':
        dat = unbin_data #for binned data for %p spec shift 
        plot_title = title + ' Flux shifted '+ str(velocity)+'km/s'
    
        #replace wavelengths with velocity shifted wavelengths 
        for i in range(len(dat['wave'])):
            dat['wave'][i] = get_obswave(dat['wave'][i], velocity)[0] #redshift
        plot_ax = plot_p_single([], [], [], dat, data, binsize, flx_adjust, plot_title)
    
    if shift_spec == 'p':
        dat = data #for binned data for %p spec shift 
        plot_title = title + '%P shifted '+ str(velocity)+'km/s'
    
        #replace wavelengths with velocity shifted wavelengths 
        for i in range(len(dat['wave'])):
            dat['wave'][i] = get_obswave(dat['wave'][i], velocity)[0] #redshift
        plot_ax = plot_p_single([], [], [], unbin_data, dat, binsize, flx_adjust, plot_title)
    
    if shift_spec == 'polflux':
        dat = data
        p = P_tot(dat['q'], dat['u'])
        perr = P_err(p, dat['q'], dat['u'], dat['qerr'], dat['uerr'])
        #pdebi = P_debi3(p, perr) #debias with flux makes eveything negative (?)
        pflux = dat['flx']*p
        plot_title = title + 'Polarized Flux'
        plot_ax = plot_p_single([], [], [], unbin_data, data, binsize, flx_adjust, plot_title)
        plot_ax.step(dat['wave'], pflux*0.003, color = 'red')
    
    """if shift_spec == 'shift_polflux': #WIP currently still calculating polfux at rest then shifting it
        #I think I would like it to shift the polarization then calculate polflux
        dat = data
        for i in range(len(dat['wave'])):
            dat['wave'][i] = get_obswave(dat['wave'][i], velocity)[0] #redshift
        p = P_tot(dat['q'], dat['u'])
        perr = P_err(p, dat['q'], dat['u'], dat['qerr'], dat['uerr'])
        #pdebi = P_debi3(p, perr) #debias with flux makes eveything negative (?)
        pflux = dat['flx']*p
        plot_title = title + 'Polarization shifted '+ str(velocity)+'km/s, Polflux Calculated'
        plot_ax = plot_p_single([], [], [], unbin_data, dat, binsize, flx_adjust, plot_title)
        plot_ax.step(dat['wave'], pflux*flx_adjust*.2, color = 'red')"""
    
    return(plot_ax)

def vel_comp_stack(unbin_data, flx_adjust, bin_data, rest_wave, velocity):
    """ Creates two figures, 1: flx, q, u stacked. 2: flx, p, pa stacked.
    Inputs are a list of unbinned data tables, list of flx adjusment factors, list of binned data tables 
    list of rest wavelengths and velocity to determine range to plot data over.
    To compare a single line over multiple epochs provide one value in rest_wave list 
    and all epoch of data in data tables and flx adjust lists. 
    To compare multiple lines within a single epoch provide all lines rest wavelengths in 
    rest_wave list and single epochs in data tables lists, in this case
    flx_adjust can be provided or set = [1] since it's all on the same scale already. 
    Examples: 
    vel_comp_stack(all_epochs_unbin, flx_adjust, all_epochs20, [6300], 20000) #single line all epochs
    vel_comp_stack(data2, [1], epoch2_20, [5876, 6678, 7065], 20000) #multiple lines one epoch"""
    
    #If comparing same line over different epochs:
    if len(rest_wave) == 1: 
        blue = get_obswave(rest_wave[0], velocity)[1] #blue observed wavelength 20000km/s out
        red = get_obswave(rest_wave[0], velocity)[0] #red observed wavelength 20000km/s out
    
        fig, ax = plt.subplots(3, 1, sharex = True, figsize = (5, 10))
        plt.subplots_adjust(wspace=0, hspace=0)
        ax[0].set_ylabel('Relative Flux')
        ax[1].set_ylabel('%q')
        ax[2].set_ylabel('%u')
        ax[2].set_xlabel('Velocity (km/s)') 
        
        
        #compare line region between epochs in velocity space
        count = 1
        for epoch, adjust in zip(unbin_data, flx_adjust):
            reg = get_lines(epoch, blue, red)
            reg_vel = get_vel_column(reg, rest_wave[0])
            ax[0].plot(reg_vel['vel'], reg_vel['flx']*adjust, label = 'epoch'+str(count))
            count = count + 1
        ax[0].legend(loc="upper right", fontsize = 12)
        
        for epoch in bin_data:
            reg = get_lines(epoch, blue, red)
            reg_vel = get_vel_column(reg, rest_wave[0])
            ax[1].step(reg_vel['vel'], reg_vel['q'])
            ax[2].step(reg_vel['vel'], reg_vel['u'])
        ax[2].set_xticks([-velocity, -velocity/2, 0, velocity/2, velocity])    

        fig2, ax2 = plt.subplots(3, 1, sharex = True, figsize = (5, 10))
        plt.subplots_adjust(wspace=0, hspace=0)
        ax2[0].set_ylabel('Relative Flux')
        ax2[1].set_ylabel('%P')
        ax2[2].set_ylabel('PA')
        ax2[2].set_xlabel('Velocity (km/s)') 
        
        #compare line region between epochs in velocity space
        count = 1
        for epoch, adjust in zip(unbin_data, flx_adjust):
            reg = get_lines(epoch, blue, red)
            reg_vel = get_vel_column(reg, rest_wave[0])
            ax2[0].plot(reg_vel['vel'], reg_vel['flx']*adjust, label = 'epoch'+str(count))
            count = count + 1
        ax2[0].legend(loc="upper right", fontsize = 12)
        
        for epoch in bin_data:
            reg = get_lines(epoch, blue, red)
            reg_vel = get_pa_column(get_vel_column(reg, rest_wave[0]))
            ax2[1].step(reg_vel['vel'], reg_vel['p'])
            ax2[2].step(reg_vel['vel'], reg_vel['pa'])
        ax2[2].set_xticks([-velocity, -velocity/2, 0, velocity/2, velocity])
        #ax2[2].set_xticks([-20000, -10000, 0, 10000, 20000])

    else:
        #compare multiple lines in same epoch
        fig, ax = plt.subplots(3, 1, sharex = True, figsize = (5, 10))
        plt.subplots_adjust(wspace=0, hspace=0)
        ax[0].set_ylabel('Relative Flux')
        ax[1].set_ylabel('%q')
        ax[2].set_ylabel('%u')
        ax[2].set_xlabel('Velocity (km/s)')
        #set rest wave 
        ax[0].axvline(x=0, c = 'k'); ax[1].axvline(x=0, c = 'k'); ax[2].axvline(x=0, c = 'k')
        fig2, ax2 = plt.subplots(3, 1, sharex = True, figsize = (5, 10))
        plt.subplots_adjust(wspace=0, hspace=0)
        ax2[0].set_ylabel('Relative Flux')
        ax2[1].set_ylabel('%P')
        ax2[2].set_ylabel('PA')
        ax2[2].set_xlabel('Velocity (km/s)')
        #set rest wave 
        ax2[0].axvline(x=0, c = 'k'); ax2[1].axvline(x=0, c = 'k'); ax2[2].axvline(x=0, c = 'k')
        
        for line in rest_wave:
            blue = get_obswave(line, velocity)[1]
            red = get_obswave(line, velocity)[0]
            reg_vel = get_pa_column(get_vel_column(get_lines(bin_data, blue, red), line)) 
            reg_flx = get_vel_column(get_lines(unbin_data, blue, red), line)
            #flx, q, u
            ax[0].plot(reg_flx['vel'], reg_flx['flx']*flx_adjust[0], label = str(line))
            ax[1].step(reg_vel['vel'], reg_vel['q'])
            ax[2].step(reg_vel['vel'], reg_vel['u'])
            #flx, p, pa
            ax2[0].plot(reg_flx['vel'], reg_flx['flx']*flx_adjust[0], label = str(line))
            ax2[1].step(reg_vel['vel'], reg_vel['p'])
            ax2[2].step(reg_vel['vel'], reg_vel['pa'])
        ax[2].set_xticks([-velocity, -velocity/2, 0, velocity/2, velocity])
        ax2[2].set_xticks([-velocity, -velocity/2, 0, velocity/2, velocity])
        ax[0].legend(loc= 'upper right'); ax2[0].legend(loc= 'upper right')
    return()

def plot_comp_Ps(data_t, title, avePA = 'none'):
    from get_data import P_err, P_debi, P_debi2, QRSP, URSP, PA, P_tot 
    
    PAs = PA(data_t['q'], data_t['u'])
    if avePA == 'none':
        avePA = np.mean(PAs) #calculate average PA from data unless a value is given 
    
    tradp = P_tot(data_t['q'], data_t['u'])
    debip = P_debi(data_t['q'], data_t['u'], data_t['qerr'], data_t['uerr'])
    QRSP = QRSP(data_t['q'], data_t['u'], avePA)
    URSP = URSP(data_t['q'], data_t['u'], avePA)
    PRSP = P_tot(QRSP, URSP)
    trad_perr = P_err(tradp, data_t['q'], data_t['u'], data_t['qerr'], data_t['uerr'])
    debi_perr = P_err(debip, data_t['q'], data_t['u'], data_t['qerr'], data_t['uerr'])
    PRSP_perr = P_err(PRSP, data_t['q'], data_t['u'], data_t['qerr'], data_t['uerr'])
    debip2 = P_debi2(tradp, trad_perr)

    fig, ax = plt.subplots(4, 1, figsize = (10, 12), sharex = True)
    fig.suptitle(title, fontsize=16, y=.92)
    plt.subplots_adjust(hspace = 0)
    ax[0].step(data_t['wave'], debip2, c = 'lime', linewidth = 8, label = 'sqrt(|p**2 - σ**2|)')
    ax[0].step(data_t['wave'], tradp, c = 'b', linewidth = 6, label = 'trad_p')
    ax[0].step(data_t['wave'], debip, c = 'y', linewidth = 4, label = 'debi_p')
    ax[0].step(data_t['wave'], PRSP, c = 'm', linewidth = 2, label = 'PRSP')
    ax[0].set_ylabel('%P')
    ax[0].legend()

    ax[1].step(data_t['wave'], debip2/trad_perr, linewidth = 8, c = 'lime')
    ax[1].step(data_t['wave'], tradp/trad_perr, linewidth = 6, c = 'b')
    ax[1].step(data_t['wave'], debip/debi_perr, linewidth = 4, c = 'y')
    ax[1].step(data_t['wave'], PRSP/PRSP_perr, linewidth = 2, c = 'm')
    ax[1].axhline(np.mean(debip2/trad_perr), c = 'purple', label = 'debi2. Ave P/σ'+str(np.round(np.mean(debip2/trad_perr), 2)))
    ax[1].axhline(np.mean(tradp/trad_perr), c = 'k', label = 'trad. Ave P/σ: ' +str(np.round(np.mean(tradp/trad_perr), 2)))
    ax[1].axhline(np.mean(debip/debi_perr), c = 'r', label = 'debi. Ave P/σ: ' +str(np.round(np.mean(debip/debi_perr), 2)))
    ax[1].axhline(np.mean(PRSP/PRSP_perr), c = 'g', label = 'PRSP Ave P/σ: ' +str(np.round(np.mean(PRSP/PRSP_perr), 2)))
    ax[1].legend()
    #ax[1].set_ylim(0, 1)
    ax[1].set_ylabel('P/σ')
    
    ax[2].step(data_t['wave'], QRSP, c = 'orange', label = 'q_RSP')
    ax[2].step(data_t['wave'], URSP, c = 'c', label = 'u_RSP')
    ax[2].legend()
    ax[2].set_ylabel('RSP')    
    
    ax[3].step(data_t['wave'], PAs, c = 'k')
    ax[3].axhline(avePA, c = 'r', label = 'Ave PA: ' +str(np.round(avePA, 2)))
    ax[3].set_ylabel("θ (°)")
    ax[3].set_xlabel("Wavelength (Å)")
    ax[3].legend()
    
    return()


def plot_pannels(unbinned_t, bin_data_t, unbin_flux_factor, bin_flux_factor, adjust_pa = 'on', title = 'none'):
    
    """Creates the classic five plots (total polarization,polarized flux, %U, %Q, PA vs wavelength) with the flux 
    plotted in a shaded region of each pannel. 
    
    Takes in two astropy tables for unbinned and binned data
    unbinned_file must contain a column call 'flx' (so make sure to pass that to get_fits if that's used)
    
    Flux factors are what to multiply the unbinned and binned flux by to be on the same scale as the polarization. 
    adjust_pa is automatically set to add 180 for all values, if you don't want this adjust_pa = 'off'.
    
    Provide a title to save plot otherwise the plot will be titled 'none' and will not be saved. 
    """
    
    import numpy as np
    import matplotlib.pyplot as plt
    
    Q = np.array(bin_data_t['q']) #Q and U should already be in percent 
    U = np.array(bin_data_t['u'])
    
    F = np.array(unbinned_t['flx'])*unbin_flux_factor #Adjust scale of relative flux to magnitude of polarizations
    
    P_Flux = bin_data_t['p']*bin_data_t['flx']*bin_flux_factor
    
    PA = np.degrees(.5*np.arctan2(U, Q)) #calculate PA
    
    if adjust_pa != 'on': #decide to add 180 to PA or not
        pass
    else:
        for i in range(len(PA)):
            PA[i] = PA[i]+180

    PA_ave = [] #creat list for plotting line of average PA 
    for i in range(len(PA)): PA_ave.append(np.nanmean(PA)) #calculate average PA
    
    fig, (ax0, ax1, ax2, ax3, ax4) = plt.subplots(5, sharex = True, figsize=(20,14))
    fig.suptitle(title, fontsize=24, y=.95)
    ax0.step(bin_data_t['wave'], bin_data_t['p'], 'r'); ax0.fill_between(unbinned_t['wave'], F, color = 'skyblue', label = 'Relative Flux'); ax0.set_ylabel("%Polarization", fontsize = 'medium'); ax0.legend(fontsize='medium')
    ax1.step(bin_data_t['wave'], P_Flux, 'r'); ax1.fill_between(unbinned_t['wave'], F, color = 'skyblue'); ax1.set_ylabel("Polarized Flux", fontsize = 'medium')
    ax2.step(bin_data_t['wave'], Q, 'r'); ax2.fill_between(unbinned_t['wave'], F, color = 'skyblue'); ax2.set_ylabel("% q", fontsize = 'medium')   
    ax3.step(bin_data_t['wave'], U, 'r'); ax3.fill_between(unbinned_t['wave'], F, color = 'skyblue' ); ax3.set_ylabel("% u", fontsize = 'medium')   
    ax4.step(bin_data_t['wave'], PA, 'r'); ax4.plot(bin_data_t['wave'], PA_ave, 'k', label = 'Average PA: '+ str(round(np.nanmean(PA))) + u'\xb0', linewidth = 3);ax4.set_ylabel("Position Angle", fontsize = 'medium'); ax4.legend(fontsize = "medium"); ax4.set_xlabel("Wavelength(Å)")
    
    plt.subplots_adjust(wspace=0, hspace=0)
    #if title != 'none':
        #fig.savefig("Pannel_Plots/"+ title)
        #print("Plot has been saved under filename " + str(title) + " in the Pannel_Plots directory")
       
    return (print("provide title to save plot"))


def best_fit(Q, U, Qsig, Usig, qmin = 'none'):
    #(X, Y, Xsig, Ysig)
    from scipy import odr 
    from scipy.optimize import curve_fit
    from sklearn.linear_model import LinearRegression
    
    """Creates an error-weighted linear best fit line for Q/U plots given 4 lists of percent. 
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
    #odr.set_iprint(init=2, iter=2, final=2) #uncomment to see what odr.run does
    out = odr.run()
    
    #fit peramaters
    popt = out.beta
    perr = out.sd_beta
    #chi_square = out.sum_square/out.iwork[10]#same as out.res_var
    chi_square  = out.res_var #not as good an estimate as using lin_bfit
    #print("fit parameter 1-sigma error")
    #print("Slope = "+ str(popt[0])+" +- "+str(perr[0]))
    #print("Intercept = " + str(popt[1]) + " +- " + str(perr[1]))
    #print("———————————–")
    
    if qmin == 'none':
        x_fit = np.linspace(min(Q), max(Q)) #x pts
    else:
        x_fit = np.linspace(qmin[0],qmin[1]) #x pts -extending the line 
    fit = func(popt, x_fit) #y pts
    
    #for 5-sigma confidence level curve:
    nstd = 5. # to draw 5-sigma intervals
    popt_up = popt + nstd * perr
    popt_dw = popt - nstd * perr
    fit_up = func(popt_up, x_fit)
    fit_dw= func(popt_dw, x_fit)
    
    #links for weighted fit:
    #https://micropore.wordpress.com/2017/02/07/python-fit-with-error-on-both-axis/
    #https://docs.scipy.org/doc/scipy/reference/odr.html
    
    return (x_fit, fit, fit_up, fit_dw, popt[0], perr[0], chi_square)

def plot_QU(binned_data, title, size = 3.5, save = 'none', epoch_labels = 'none'):
    
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
    ax.axis('square')
    rcParams["font.size"]= 20
    fig.suptitle(title, fontsize=20, y=.85, x=.44) 
    ax.axhline(0, color = 'k', linewidth = 1, zorder = 1)
    ax.axvline(0, color = 'k', linewidth=1, zorder = 2)
    ax.plot(Q, U, c = 'lightgrey', linewidth = 1, zorder = 3 )
    #ax.axis('square'); ax.set_xlim([-float(size), float(size)]); ax.set_ylim([-float(size), float(size)]) #comment out to see zoomed in 
    ax.set_xlim(-3, 2.5); ax.set_ylim(-3, 4.5); #ax.axis('square')
    ax.set_xlabel("% q Polarization", fontsize=18); ax.set_ylabel("% u Polarization", fontsize=18)

     
    if epoch_labels == 'none':
        #if plotting a a single epoch of data plot colors accoring to wavelength
        ax.errorbar(Q,U, xerr = Qsig, yerr = Usig) #plot individual point uncertainties 
        im = ax.scatter(Q,U, marker = 'o', c=wl, s=20, zorder = 4, cmap = 'turbo') #c=color set to a different color for each point in wavelength array, s=size 
        fig.colorbar(im).ax.set_ylabel('Wavelength(Å)', fontsize=20) #shows colorbar and labels it
        ax.plot(x_fit, fit, "k", lw=2, label = 'best fit', zorder = 5) #best fit line
        #ax.errorbar((-float(size) + .005), (float(size) - .0058), xerr = Qave_err, yerr = Uave_err) #add error bar example in location based on grid size 
        ax.add_patch(plt.Circle((0, 0), radius = .66, edgecolor='red', facecolor = 'none', label = 'ISP'))
        #ax.text((-float(size)+ .008), (float(size) - .0055), "Average Error", fontsize = 13)
        ax.errorbar(2,-.4, xerr = Qave_err, yerr = Uave_err, c ='purple')
        ax.text(2.1,-.4, "Ave Error", fontsize = 13)
        
    else:
        #if ploting average data points (wavelength doesn't matter) plot all points as same color 
        im = ax.scatter(Q,U, marker = 'o', s=20, zorder = 5)
        ax.errorbar(Q, U, yerr=Usig, xerr=Qsig, ecolor="grey", hold=True, fmt="none", zorder = 4 ) #comment out to see points by wavelength
        for i, txt in enumerate(epoch_labels):
            ax.annotate(txt, xy = (Q[i]+0.0005, U[i]), fontsize = 14) #label each data point according to epoch_labels provided     
    
    ax.legend(loc="upper right",fontsize='x-small')
    if save != 'none':
        fig.savefig("QU_Plots/"+ title)
        print("figure saved in QU_Plots directory under" + str(title))
    else:
        print("set save = on to save figure")
    fig.show()
     
    #return (print("Q vs U plot has been saved under filename " + str(title) + " in the QU_Plots directory"))
    return() 
   
def plot_line_velspace(data_t, rest_wave, velocity, flx_factor, epoch, save = "off", plot = "off", region = "full"):
    """Creates a velocity space plot centered on desired line. 
    Takes in a data table (from Bin_data or get_txtFITS), the wavelength of the line to center plot on, 
    the velocity range to go out to on either end, the flux factor so flux fits to scale with polarization, 
    an epoch for the plot title (no spaces) and option to save if save = on to folder VelSpace_plots/epoch+rest_wave
    (save is automatically off so plot will not save)
    if region = 'blue' only the blue side of the line will be ploted and returned, default is "full" so red and blue centered on the line"""

    title = str(epoch) + " Velocity Region Centered on " + str(rest_wave)
    vel=[] #creat list for velocity values
    wave_range = get_obswave(rest_wave, velocity) #calculate wavelengths on either side of rest that pertain to velocity region desired
    if region == 'blue':
        line_data = get_lines(data_t, wave_range[1], rest_wave) #get only blue end of line region
    else:
        line_data = get_lines(data_t, wave_range[1], wave_range[0]) #get lines from data table in between wavelengths in velocity range
    #print(line_data)
    line_data = get_vel_column(line_data, rest_wave)
    vel = line_data['vel']
    
    if plot != "off":
        #plot total polarization and flux in velocity range
        plt.figure(figsize= (10,7)); plt.title(title, fontsize = '20')
        plt.step(vel, line_data['p'], color = 'k', label = "Total Polarization")
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

def plot_all_flx(unbin_epoch_tables_list, flx_adjust_list, epoch_names):
    """ Creates a figure of stacked flux for multiple epochs given a list of unbinned data tables, 
     a list of factors to multiply flx by to get them all on the same scale and a list of epoch names.
     Returns a stacked flux on an axis at which point individual features on plot can be labeled from"""
    
    # loop through the length of list of unbinned epoch data tables and keep track of index
    fig, (ax) = plt.subplots(1, sharex = True, figsize = (14, 10), dpi = 1200)
    ax.set_xlabel('Rest Wavelength (Å)', fontsize = 20)
    #ax.set_xticks(np.arange(3500, 8400, 500))#show only x-axis for last plot
    ax.set_ylabel('Relative Flux', fontsize = 20)
    ax.tick_params(axis='x', labelsize=20); ax.tick_params(axis='y', labelsize=20)
    
    #count number of epochs
    num_epochs = 0
    for count, item in enumerate(unbin_epoch_tables_list):
        num_epochs = num_epochs +1
    space = num_epochs*4
    
    for n, epoch in enumerate(unbin_epoch_tables_list):
        count = n

        ax.plot(epoch['wave'], epoch['flx']*flx_adjust_list[count]+space, color = 'k')#plot flux with adjusted scale and inline with polarization
        #ax.plot(epoch['wave'], epoch['flx']*+space, color = 'k')#plot flux with adjusted scale and inline with polarization
        
        #ax.annotate(epoch_names[count], (epoch['wave'][-1] - 300, epoch['flx'][-1]*flx_adjust_list[count]+space+1), fontsize = 16)
        space = space - 2
  
    return(ax)
def plot_all_epoch_pannels(unbin_data_list, binned_data_list, flx_adjust_list, epoch_names, plot_data = 'p', line_list = 'none', color_list = 'none', label_list = 'none', regions_list = 'none', step_color = '#984EA3'):
    """ Creates stacked pannels of flx and either %p or PA for each epoch.
    Inputs required are: 
    a list of unbinnned data tables, a list of binned data tables, a flux adjust list, a list of epoch names
    Options are:
    If PA is the desired second plot line specify plot_data = 'pa' and make sure the flx_adjust list  is made up of value to match the scale of the angles.
    To plot single vertical lines marking features provide a list of lines, colors for those lines and labels.
    To shade in vertical region for continuum estimate provide list of regions_list = [[region1_b_wl, region1_r_wl], [region2_b_wl, region2_r_wl]...]
    Can specify step_color to change color the polarization spectra is plotted in. 
    The first item returned is the ax, the second is a table of averages for the continuum regions provided by epoch with regions stacked on eachother. 
        ie: output = plot_all_epoch_pannels(blah, blah, regions_list = [[5000, 5100], [7400, 7500]]) -> output[1] = Table[epoch1_reg1aves, epoch2_reg1aves, epoch1_reg2_aves, epoch2_reg2aves]
    """
    num_epochs = 0
    for count, item in enumerate(unbin_data_list):
        num_epochs = num_epochs +1
    
    fig, ax = plt.subplots(num_epochs, 1, sharex = True,  figsize = (15, 12), dpi = 1200)
    fig.subplots_adjust(right = .85,  wspace=0, hspace=0)
    
    # loop through the length of list of unbinned epoch data tables and keep track of index
    for n, epoch in enumerate(unbin_data_list):
        # add a new subplot iteratively
        ax[n] = plt.subplot(num_epochs, 1, n+1)#stack them
        #ax[n].set_xlim(3500, 7900); #
        ax[n].set_xticks(np.arange(4000, 8000, 250)); ax[n].set_xticklabels([]) #share same x-axis scale but don't show for each of them 
        ax[n].xaxis.set_major_locator(MultipleLocator(250)); ax[n].xaxis.set_minor_locator(MultipleLocator(50))
        
        ax[n].plot(epoch['wave'], epoch['flx']*flx_adjust_list[n], c = 'k') #label each one by epoch
        ax[n].tick_params(axis='x', which='major', labelsize=14)
    #plot epoch labels
        if plot_data == 'pa':
            ax[n].annotate(epoch_names[n], (4200, 2.8), fontsize = 16)
        else:
            ax[n].set_ylim(0.001, 3.7);
            ax[n].yaxis.set_major_locator(MultipleLocator(1)); ax[n].yaxis.set_minor_locator(MultipleLocator(.2))
            ax[n].tick_params(axis='both', which = 'both', labelleft = True, labelsize=14, top = True, right = True)
    
            if n >= 4:
                ax[n].annotate(epoch_names[n], (4050, 4.2), fontsize = 16)
                ax[n].set_ylim(0, 5.8); ax[n].yaxis.set_major_locator(MultipleLocator(2)); ax[n].yaxis.set_minor_locator(MultipleLocator(.2))
                ax[n].tick_params(axis='y', which = 'both', labelleft = True, labelsize=14)
            else:
                ax[n].annotate(epoch_names[n], (4050, 3), fontsize = 16)
        
        #plot spectra feature
        if line_list != 'none':
            #if regions == 'off':
            #plot single lines, colorcoded and labeled in top pannel
            #if n == 5:     #uncomment if n == 5: statement for 2012au (used to skip epoch 6)
                #continue
            for xc, c in zip(line_list, color_list):
                ax[n].axvline(x=xc, c=c, alpha = 0.5) #plot veritcle line at location xc in linelist with color c in colorlist
            if n == 0:
                for xc, label in zip(line_list, label_list):
                    ax[n].text(xc + 20, 3.1, label, fontsize = 16) #add label at verticle position pmax next to line at xc, with label lab from label list 
              
    #highlight continuum regions (if lists of region bounds are provided)
    ave_tab = Table() #return table of averages calculated between regions
    if regions_list != 'none':
        count_lists = 0 #index first list in list of lists
        for region in regions_list:
            #print(region)
            aves = get_aves(binned_data_list, region) #calculate averages for data in region
            #print(aves)
            ave_tab = vstack([ave_tab, aves])
            for n, epoch in enumerate(binned_data_list):
                #if n == 5: #    uncomment if n == 5: statement for 2012au (used to skip epoch 6)
                   # continue
                ax[n].axvspan(region[0], region[1], color = 'lightgrey') #highlight region 
            #ax[n].annotate(ave_p + "% \n \u00B1" + ave_perr, xy =(region[0] + 50, space+8 ), fontsize = 9) #lable polarization and error in region where xy = is a really complicated way of locating the label        
        count_lists = count_lists + 1 #move to next region list index    
     
    for n, epoch in enumerate(binned_data_list):#plot polarization spectra and errors
        # add a new subplot iteratively
        ax[n] = plt.subplot(num_epochs, 1, n+1)#stack them
        ax[n].set_xlim(4000, 8000); #ax[n].set_xticks([])#share same x-axis scale but don't show for each of them
        
        if plot_data == 'p':
            ax[n].step(epoch['wave'], epoch['p'], c = step_color)
            #ax[n].fill_between(epoch['wave'], epoch['p']-epoch['perr'],  epoch['p']+epoch['perr'], color = 'lightgrey', zorder = 1)

            #ax.axhline(y = 0.67, c = 'r')
            if n == len(epoch_names)/2-1:   
                ax[n].set_ylabel('Percent Polarization', fontsize = 16)
            for x, val in enumerate(epoch['perr']):
                if val < 3:
                    ax[n].errorbar(epoch['wave'][x], epoch['p'][x], yerr = val ,c = 'thistle', elinewidth = 4, capsize=2, zorder = 1)
            
        if plot_data == 'pa':
            pa = PA(epoch['q'], epoch['u'], adjust_pa = 'off')
            ax[n].step(epoch['wave'], pa)
            if n == len(epoch_names)/2-1:   
                ax[n].set_ylabel('Position Angle (°)', fontsize = 12)
                
        if plot_data == 'pflux':
            ax[n].step(epoch['wave'], epoch['p']*epoch['flx']*flx_adjust_list[n])
            if n == len(epoch_names)/2-1:   
                ax[n].set_ylabel('Polarized Flux', fontsize = 12)
        
        if plot_data == 'q':
            ax[n].set_ylim(-3.7, 3.7)
            ax[n].yaxis.set_major_locator(MultipleLocator(1)); ax[n].yaxis.set_minor_locator(MultipleLocator(.2))
            ax[n].tick_params(axis='both', which = 'both', labelleft = True, labelsize=14, top = True, right = True)   
            ax[n].step(epoch['wave'], epoch['q'], c= '#377EB8')
            if n == len(epoch_names)/2-1:   
                ax[n].set_ylabel('q(%)', fontsize = 12)
            for x, val in enumerate(epoch['qerr']):
                if val < 3:
                    ax[n].errorbar(epoch['wave'][x], epoch['q'][x], yerr = val ,c = 'thistle', elinewidth = 4, capsize=2, zorder = 1)
        
        if plot_data == 'u':
            ax[n].set_ylim(-3.7, 3.7)
            ax[n].yaxis.set_major_locator(MultipleLocator(1)); ax[n].yaxis.set_minor_locator(MultipleLocator(.2))
            ax[n].tick_params(axis='both', which = 'both', labelleft = True, labelsize=14, top = True, right = True)
            ax[n].step(epoch['wave'], epoch['u'], c = '#4daf4a')
            if n == len(epoch_names)/2-1:   
                ax[n].set_ylabel('u(%)', fontsize = 12)
            for x, val in enumerate(epoch['uerr']):
                if val < 3:
                    ax[n].errorbar(epoch['wave'][x], epoch['u'][x], yerr = val ,c = 'thistle', elinewidth = 4, capsize=2, zorder = 1)   
            
    #show only x-axis for last plot
    ax[n].set_xlabel('Rest Wavelength  (Å)', fontsize = 16)
    ax[n].set_xticks(np.arange(4000, 8000, 250)); ax[n].set_xticklabels(np.arange(4000, 8000, 250))
    ax[n].xaxis.set_major_locator(MultipleLocator(250)); ax[n].xaxis.set_minor_locator(MultipleLocator(50))
    return(ax, ave_tab)

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
        image.axvline(x=...)"""
  
    
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

def QU(binned_data, ax = None, bfit = 'on', size = 4.3, color_data = 'wl', epoch_labels = 'none', cmap_choice = 'turbo', bfit_c = '#999999', regions = 'none', marker = 'o', ec = 'face'):
    """Returns QU plot given a binned_data table. The return object should be added to a figure in order
    to do fig.colorbar(QU_return) and fig.suptitle to make it look pretty. Multiple QU plots can be added 
    as subplots of a figure by: 
        import matplotlib.pyplot as plt
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize = (12, 5))
        e1 = QU(epoch1_20, 'E1', ax1)
        e2 = QU(epoch2_20, 'E2', ax2)
        fig.colorbar(e1).set_label('Wavelength(Å)') 
    A best fit line can be calculated outside the function and passed to it to exclude points that are in the QU plot. 
    If bfit is not specified the bestfit in calculated in the function using all the QU points. 
    If a list of values specifying the blue end, rest wave and red end of a feature are specified in in regions = [] then 
    the points are colorcoded accordingly.
    If color_data is given as a list of ranges this list will be used to colorcode points:
        ex: if say multiple epochs have different wavelength ranges but are going to be plotted as subplots together they 
        can share a common color bar if the same list spanning all the ranges of data is given for color_data"""
    
    if ax == None:
        ax = plt.gca()
        ax.set_xlim([-float(size), float(size)+0.1]); ax.set_ylim([-float(size), float(size)]) #comment out to see zoomed in  
    
    ax.axhline(0, color = 'k', linewidth = 1, zorder = 1)
    ax.axvline(0, color = 'k', linewidth=1, zorder = 2)
   #ax.set_xlim([-size, size]); ax.set_ylim([-size, size])
    #ax.axis('square'); 
    ax.set_aspect('equal', 'box')
    
    wl = binned_data['wave']
    #convert fraction to percentage and ignore nan values
    Q = [x for x in binned_data['q'] if np.isnan(x) == False] 
    U = [x for x in binned_data['u'] if np.isnan(x) == False]
    Qsig = [x for x in binned_data['qerr'] if np.isnan(x) == False]
    Usig = [x for x in binned_data['uerr'] if np.isnan(x) == False]
    
    wl, Q = zip(*zip(wl, Q)) #makes sure wavelength and other lists are the same length for plotting
                             #by pairing values up in each then dropping unpaired values 
    
    if bfit != 'off':
        if bfit == 'on':
            bfit = best_fit(Q, U, Qsig, Usig, qmin=[-size, size])
            x_fit = bfit[0]
            fit = bfit[1]
        else:
            x_fit = bfit[0]
            fit = bfit[1]
        
        ax.plot(x_fit, fit, bfit_c, lw=2, zorder = 5) #best fit line
        #ax.set_xlim([-size, size]); ax.set_ylim([-size, size])
    
    #calculate average error
    Qave_err = np.nanmean(Qsig)
    Uave_err = np.nanmean(Usig)
    #print("AVE ERRORS: " + str(Qave_err) +  str(Uave_err))
    if Qave_err < 0.05:
        Qave_err = 0.05
    if Uave_err < 0.05:
        Uave_err = 0.05 

    
    if ec != 'face':
    #if marker != 'o':
        #ax.errorbar((-float(size) + .5), (float(size) - .5), xerr = Qave_err, yerr = Uave_err, color = ec) #add error bar example in location based on grid size 
        ax.plot(Q, U, c = ec, linewidth = 2, zorder = 3 )
        im = ax.scatter(Q,U, marker = marker, c = color_data, edgecolors=ec, s=75, zorder = 4, cmap = cmap_choice, label = ' ')#, norm=DivergingNorm(0)) #for spec_comp to scale evenly
        #ax.errorbar(Q,U, xerr = Qsig, yerr = Usig)
        
    else :
        ax.errorbar((-float(size) + .5), (float(size) - .5), xerr = Qave_err, yerr = Uave_err, color = 'k') #add error bar example in location based on grid size 
        im =ax.plot(Q, U, c = '#999999', linewidth = 1, zorder = 3 ) #'#999999' = grey

    #ax.set_xlabel("q(%)", fontsize = 18); ax.set_ylabel("u(%)", fontsize = 18); ax.yaxis.set_label_coords(-.13, .5)
    #ax.tick_params(labelsize=18)
    
    if epoch_labels == 'none':
        #if plotting a a single epoch of data plot colors accoring to wavelength
        #if another column of data is provided for colorcoding use that (for example 'vel' column from get_vel_column to colorcode in velocity space)
        if color_data != 'wl':
            if color_data == 'constant':
                #ax.errorbar(Q,U, xerr = Qsig, yerr = Usig)
                im = ax.scatter(Q,U, marker = marker, c=bfit_c, s=20)
            if regions != 'none':
                #cmap, norm = mcolors.from_levels_and_colors(regions, ['blue', 'red'])
                im = ax.scatter(Q,U, marker = marker, s=20, zorder = 4, c = color_data, cmap = cmap_choice)#norm = DivergingNorm(0)) #c=color set to a different color for each point in wavelength array, s=size 
            #else:
                #for color by wavelengths
                #im = ax.scatter(Q,U, marker = marker, c=binned_data['wave'], s=20, zorder = 4, cmap = cmap_choice, vmax = max(color_data), vmin = min(color_data)) #if color_data is specified by a list of ranges use the min and max of that list to define colorcoding of points (so that different length data sets can share the same colorbar)
                #for colors by velocities 
                #im = ax.scatter(Q,U, marker = marker, c=binned_data['vel'], s=20, zorder = 4, cmap = cmap_choice, vmax = max(color_data), vmin = min(color_data)) #if color_data is specified by a list of ranges use the min and max of that list to define colorcoding of points (so that different length data sets can share the same colorbar)
        else:
            #ax.errorbar(Q,U, xerr = Qsig, yerr = Usig)
            im = ax.scatter(Q,U, marker = marker, c=wl, s=20, zorder = 4, cmap = cmap_choice)# norm=DivergingNorm(0)) #c=color set to a different color for each point in wavelength array, s=size 
        #fig.colorbar(im).ax.set_ylabel('Wavelength(Å)', fontsize=12) #shows colorbar and labels it
        #ax.text((-float(size)+ .8), (float(size) - .55), "Average Error", fontsize = 10)
        #ax.legend(loc="lower left", fontsize='medium')
    else:
        
        #if ploting average data points (wavelength doesn't matter) plot all points as same color 
        im = ax.scatter(Q,U, marker = marker, s=20, zorder = 5)
        #ax.errorbar(Q, U, yerr=Usig, xerr=Qsig, ecolor="grey", fmt="none", zorder = 4 ) #comment out to see points by wavelength
        for i, txt in enumerate(epoch_labels):
            ax.annotate(txt, xy = (Q[i]+0.0005, U[i]), fontsize = 14) #label each data point according to epoch_labels provided     
    
    return(im)

def plot_QU_spec_comp2(unbindata, bindata, flx_adjust, region_lists, region_names, comp_color = 'whole_region', rest_qu = 'none', blue_2reg = 'none'):
    """Creates a plot with the flux and polarization spectrum of a single epoch of data as the bottom pannel
    and two square plots comparing Q/U data of selected 2 regions in that same epoch above (hence comp2). 
    Q/U plots are color coded in velocity space with the colorbar dictated by the second region given. 
    Inputs required are: unbinned data table from get_fits, bindata table from Bin_data, 
    flx_adjust = a single number for flx adjust to adjust flx spectrum to polarization level for bottom pannel, 
    region_lists = a list of two lists for the two line regions to compare where in each list the entries are 
    [blue_limit_wave, rest_wave, red_limit_wave], region_names = a list of region names as strings and the region, 
    comp_color = a string value (blue or red) that specifies what side of the rest wavelength you want to compare in Q/U space.
    If no color is given the default is to compare the whole line region from the blue limit to the red limit
    If a list is provided for rest_qu the q and u values for the rest wavlength will be plotted as a black X on the top two panels:
        [rest_q_value_region1, rest_u_value_region1, rest_q_value_region2, rest_u_value_region2]
      
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
    ax1.set_ylabel('Percent Total Polarization', fontsize = 14)
    ax1.set_xlabel('Wavelength (Å)', fontsize = 14)
    ax1.tick_params(axis='both', which='major', labelsize=14)
    
    #top sqaures - line region QUs
    ax2.set_ylabel("% u", fontsize = 14)
    ax2.set_xlabel("% q", fontsize = 14)
    ax2.xaxis.set_label_coords(.5, -.1)
    ax2.add_patch(plt.Circle((0, 0), radius = .34, edgecolor='red', facecolor = 'none'))#isp
    ax2.tick_params(axis='both', which='major', labelsize=14)
    ax2.annotate(region_names[0], (1.2, 3.5))
    

    #ax3.set_ylabel("% u", fontsize = 14)
    ax3.set_xlabel("% q", fontsize = 14)
    ax3.xaxis.set_label_coords(.5, -.1)
    ax3.add_patch(plt.Circle((0, 0), radius = .34, edgecolor='red', facecolor = 'none'))#isp
    ax3.tick_params(axis='both', which='major', labelsize=14)
    ax3.set_yticks([])#share u axis
    ax3.annotate(region_names[1], (1.2, 3.5))
    
    
    #plot rest wavelength values of Q/U for each region
    if rest_qu != 'none':
        ax2.plot(rest_qu[0], rest_qu[1], marker = 'X', markersize = 12, color = 'k', zorder = 4)
        ax3.plot(rest_qu[2], rest_qu[3], marker = 'X', markersize = 12, color = 'k', zorder = 4)
                 
    count = 0
    for region, ax in zip(region_lists, ax_list):

        ax1.axvspan(region[0], region[1], color = 'paleturquoise', alpha = .5)
        ax1.axvspan(region[1], region[2], color = 'lightcoral', alpha = .4)
        
        #get QU data in velocity space for blue region of lines
        if comp_color == 'blue':
            if blue_2reg != 'none':
                #superblue region 
                reg1 = get_lines(bindata, blue_2reg[count][0], blue_2reg[count][1])
                reg1_vel = get_vel_column(reg1, blue_2reg[count][2]) #calculate velocity w.r.t. rest 
                reg1_qu = QU(reg1_vel, ax = ax, color_data = reg1_vel['vel'], cmap_choice = 'winter')
                if count == 0:
                    cbar_ax = fig.add_axes([0, 0.5, 0.4, 0.05])
                    fig.colorbar(reg1_qu, cax=cbar_ax, orientation='horizontal') 
                
                #blue region 
                reg2 = get_lines(bindata, blue_2reg[count][1]-20, blue_2reg[count][2])
                reg2_vel = get_vel_column(reg2, blue_2reg[count][2])
                reg2_qu = QU(reg2_vel, ax = ax, color_data = reg2_vel['vel'], cmap_choice = 'spring')
                if count == 0: #share colorbar based on first line region 
                    plt.subplots_adjust(left=0, bottom=0, right=.85, top=.99, wspace=0.42, hspace=0.43)
                    cbar_ax = fig.add_axes([.85, 0.5, 0.02, 0.5])
                    fig.colorbar(reg_qu, cax=cbar_ax).set_label('Velocity(km/s)', rotation= 90, fontsize = 14)
                
            else: #when comparing the blue side of data make the color_map green/blue
                if count == 0: 
                    reg = get_lines(bindata, region[0], region[1])
                    reg_vel = get_vel_column(reg, region[1])
                    reg_qu = QU(reg_vel, ax = ax, color_data = reg_vel['vel'], cmap_choice = 'winter')
                    #cbar_ax = fig.add_axes([.5, 0.5, 0.02, 0.5])
                    #fig.colorbar(reg_qu, cax=cbar_ax).set_label('Velocity(km/s)', rotation= 90, fontsize = 14)
                else: 
                    reg = get_lines(bindata, region[0], region[1])
                    reg_vel = get_vel_column(reg, region[1])
                    reg1_color= get_lines(bindata, region_lists[0][0], region_lists[0][1]) #calculate color range based on velocities of first region given 
                    reg1_vel_color = get_vel_column(reg1_color, region_lists[0][1])
                    if np.min(reg_vel['vel']) < np.min(reg1_vel_color['vel']): # use first region to color code incase the second region is too far red that the whole velocity region doesn't have data 
                        reg_qu = QU(reg_vel, ax = ax, color_data = reg1_vel_color['vel'], cmap_choice = 'winter')
                        cbar_ax = fig.add_axes([.85, 0.5, 0.02, 0.5])
                        fig.colorbar(reg_qu, cax=cbar_ax).set_label('Velocity(km/s)', rotation= 90, fontsize = 14)
                    else: 
                        reg_qu = QU(reg_vel, ax = ax, color_data = reg_vel['vel'], cmap_choice = 'winter')
                        cbar_ax = fig.add_axes([.85, 0.5, 0.02, 0.5])
                        fig.colorbar(reg_qu, cax=cbar_ax).set_label('Velocity(km/s)', rotation= 90, fontsize = 14)
                    
                
            
        #get QU data in velocity space for red region of lines
        if comp_color == 'red':
            reg = get_lines(bindata, region[1], region[2])
            reg_vel = get_vel_column(reg, region[1])
            reg_qu = QU(reg_vel, ax = ax, color_data = reg_vel['vel'], cmap_choice = 'hot')
        
        if comp_color == 'whole_region':
            reg = get_lines(bindata, region[0], region[2])
            reg_vel = get_vel_column(reg, region[1])
            reg_qu= QU(reg_vel, ax = ax, color_data = reg_vel['vel'], cmap_choice = 'seismic', regions = [reg_vel['vel'][0], 0, reg_vel['vel'][-1]])
            if count == 0: #share colorbar based on first line region 
                plt.subplots_adjust(left=0, bottom=0, right=.85, top=.99, wspace=0.42, hspace=0.43)
                cbar_ax = fig.add_axes([.85, 0.5, 0.02, 0.5])
                fig.colorbar(reg_qu, cax=cbar_ax).set_label('Velocity(km/s)', rotation= 90, fontsize = 14)
            """
            if count == 0:
                cbar_ax = fig.add_axes([.35, 0.5, 0.02, 0.5])
                fig.colorbar(reg_qu, cax=cbar_ax).set_label('Velocity(km/s)', rotation= 90)
            if count == 1:
                cbar_ax = fig.add_axes([.84, 0.5, 0.02, 0.5])
                fig.colorbar(reg_qu, cax=cbar_ax).set_label('Velocity(km/s)', rotation= 90)"""
        
        count = count +1
    plt.subplots_adjust(hspace=0.2, wspace = -.2)
        
    return()

def region_info(epoch_unbin, epoch_data, flx_adjust, rest_wave, vel_range, cont_val, cont_err, binsize, title = 'none'):
    print('*****************************************************************************')
    print('Continuum level given: ' + str(cont_val))
    print('*****************************************************************************')
    #convert desired region to velocity 
    vel_waves = get_obswave(rest_wave, vel_range)
    blue_wave = vel_waves[1]
    blue_data = get_lines(epoch_data, blue_wave, rest_wave)
    #calculate integrated polarization
    add_val = []
    for val in blue_data['p']:
        if val < (cont_val - cont_err):
            print("ERROR: Polarization value bellow continuum level: " + str(val) + "+/- "+ str(cont_err))
        else:
            add_val.append(val*binsize - cont_val*binsize)        
    sum_p = np.nansum(add_val)

    print('*****************************************************************************')
    print("Integrated polarization for selected region: " + str(sum_p))
    print('*****************************************************************************')
    
    #Create 3 pannel plot (colors corrospond to the same points between wavelength, qu and velocity)
    vel_blue_data = get_vel_column(blue_data, rest_wave)
    fig, ax = plt.subplots(3, 1, gridspec_kw=dict(height_ratios=[1, 3, 1]), figsize = (8, 12))
    ax[0].set_box_aspect(1/3)
    ax[1].set_box_aspect(1)
    ax[2].set_box_aspect(1/3)
    #top pannel is the entire spectra with region shaded
    ax[0].fill_between(epoch_unbin['wave'], epoch_unbin['flx']*flx_adjust, color = 'silver')
    ax[0].step(epoch_data['wave'], epoch_data['p'], c = 'k')
    ax[0].scatter(vel_blue_data['wave']-binsize, vel_blue_data['p'], c =vel_blue_data['vel'],  cmap = 'winter', marker = "_", linewidths = 2.5, s = binsize,  zorder = 3) #colorcode region to velocity from rest
    ax[0].axvspan(blue_wave, rest_wave, color = 'paleturquoise', alpha = .5)
    ax[0].set_xlabel('Wavelength (\aa)'); ax[0].set_ylabel('%P')
    ax[0].set_ylim(-1, 7)
    #middle pannel is a square QU plot of region
    ax[1].set_xlabel('% Q'); ax[1].set_ylabel('% U')
    QU(vel_blue_data, ax = ax[1], cmap_choice= 'winter')
    #bottom pannel is region in velocity space
    ax[2].set_xlabel('Velocity (km/s)'); ax[2].set_ylabel('% P')
    ax[2].step(vel_blue_data['vel'], vel_blue_data['p'], c = 'k')
    lc_marker_shift = vel_range/(len(vel_blue_data)+binsize)
    ax[2].scatter(vel_blue_data['vel']-lc_marker_shift, vel_blue_data['p'], c =vel_blue_data['vel'],  cmap = 'winter', marker = "_", linewidths = 4, s = lc_marker_shift/2,  zorder = 3)#colorcode region by velocity from rest
    blue_unbin = get_lines(epoch_unbin, blue_wave, rest_wave)
    vel_blue_unbin = get_vel_column(blue_unbin, rest_wave)
    ax[2].fill_between(vel_blue_unbin['vel'], blue_unbin['flx']*flx_adjust, color = 'silver')
    ax[2].axhspan(0, cont_val, color = 'lightcoral', alpha = .4)
    ax[2].set_ylim(-1, np.max(vel_blue_data['p']+1))
    
    if title != 'none':
        fig.suptitle(title, y = .95, fontsize = 22)
    
    return(sum_p, ax)

def plot_polar(epoch_data, lines, colors, cont_p , contPAs, contPA_errs, ax = None, vel_range = [-18000, 0], isolate = [], binsize = 15, bellow_cont = 'on', plot_cont = 'off', adjust_pa = 'on'):
    """creates 0-180degree polar plot for lines in single epoch of data. 
    The velocities verses the PA are plotted for each line given with 
    the length notating the error in PA and the width being a single bin size converted to velocity. 
    Data with polarization values bellow the given continuum level is plotted with a red outline by default unless bellow_cont is set to 'off'. 
    Data correpsonding to the maximum polarization value is plotted a solid region instead of translucent. 
    Data corresponding to the minimum flux is plotted with a lime green outline.
    Inputs are:
    epoch_data = a single epoch data table 
    lines = [list of line rest lengths]
    colors = [list of colors to color code lines by]
    cont_p = simgle continnum value
    contPAs = [list of continuum average PAs - one for each region]
    contPA_errs = [list of continuum average PA_errs - one for each region]
    ax = an axis to plot on 
    vel_range = velocity range that to look at line in by default it is set -20000 to 0 (only blue side of rest wavelength)
    isolate = alternative to plot data by list of wavelength pair lists [[blue_end_line1, red_end_line1], [blue_end_line2, red_end_line2]] 
        by default these lists are empty and a velocity range is used to specify data range to plot
    
    ***NOTE: to create an array of subplots the polar projection must be specified as follows:
            fig, ax = plt.subplots(2, 3, figsize = (30, 16), subplot_kw=dict(polar=True))
    ***NOTE: to avoid repition labeling key is not provided but should specified on figure axis using ax.annotate"""
    
    if ax == None: #if ax is not specified it will creat axis 
        fig = plt.figure(figsize = [10,8])
        ax = fig.add_subplot(projection='polar')
    ax = ax  
    
    adjust_pa = adjust_pa #default is 'on' specifies adjust PA = PA+180, 'negative' adds 180 to negative PA's only, adjust_pa = 'off' is the raw calculated PA
    
    ####Axis for data 0-180
    #ax.set_thetamin(0)
    #ax.set_thetamax(180)
    #ax.set_xticks(np.pi/180. * np.linspace(180,  0, 5)) #specify theta labels (45, 90, 135, 180)
    #ax.set_xticklabels(("180°", "135°", "90°", "45°", "0°"), fontsize = 14)
    
    ####Axis for data 90-270 becuase 180 was added to all PA's (adjust_pa = 'on')
    ax.set_thetamin(90)
    ax.set_thetamax(270)
    ax.set_xticks(np.pi/180. * np.linspace(270,  90, 5)) #specify theta labels (45, 90, 135, 180)
    ax.set_xticklabels(("270°", "225°", "180°", "135°", "90°" ), fontsize = 14)
    #ax.set_xticklabels(("360°", "315°", "270°", "225°", "180°" ), fontsize = 14)
    
    ax.set_ylim([vel_range[1], vel_range[0]])
    ax.set_yticks(np.linspace(vel_range[0], vel_range[1], 11)); 
    #ax.set_yticklabels(("-18000", "-14000", "-10000", "-5000", "0"), fontsize = 14)
    
    #ax.set_yticklabels(np.linspace(-20000, 0, 8, endpoint=False), fontsize = 18)
    ax.set_theta_zero_location("E") #orientation of plot
    #plot continuum as a shaded slice
    if plot_cont != 'off':
        for contpa, conterr in zip(contPAs, contPA_errs):
            if adjust_pa == 'on':
                contpa = contpa + 180 #for adjust_pa = 'on' axis 90-270
            cont_minPA = np.deg2rad(contpa - conterr) #account for PA errors
            cont_maxPA = np.deg2rad(contpa + conterr)
            cont_PAs = np.linspace(cont_minPA, cont_maxPA, 100) #get all angles between PA limits
            ax.fill_between(cont_PAs,  vel_range[0], vel_range[1], color = '#999999', alpha = 0.5)


    if not(isolate):
        #if regions to plot are not specified by isolated wavelengths use a given velocity space     
        for rest_wave, color in zip(lines, colors):
            bluelim_wave = get_obswave(rest_wave, vel_range[0])[0]#get wavelength for blue end of velocity space
            redlim_wave = get_obswave(rest_wave, vel_range[1])[1]# get wavelength for red end of velocity space
            reg = get_lines(epoch_data, bluelim_wave, redlim_wave) #get data for region 
            reg = get_vel_column(reg, rest_wave) #get velocity column 
            reg = get_pa_column(reg, adjust_pa) #get pa and pa_err column ('on' specifies adjust PA = PA+180)
    
            pmax = np.nanmax(reg['p'])
            fmin = np.nanmin(reg['flx'])
            for data_line in reg: #plot data above continuum level
                if data_line['p'] >= cont_p:
                    pa = data_line['pa']; pa_error = data_line['paerr']; vel = data_line['vel']
                    theta_width = np.deg2rad((pa + pa_error) - (pa - pa_error)) #width of line is error of PA
                    height = get_velocity(rest_wave, rest_wave+binsize) #thickness of line is one binsize converted to velocity
                    if data_line['p'] == pmax: #plot data with pmax as a solid color rather than translucent with black outline
                        im = ax.bar(x=np.deg2rad(pa), height = height, width=theta_width, bottom=vel, color = color, edgecolor = 'k', linewidth = 4)
                    if data_line['flx'] == fmin: #plot data corresponding to flux min with lime green outline and translucent so it does not come out darker 
                        im = ax.bar(x=np.deg2rad(pa), height = height, width=theta_width, bottom=vel, color = color, edgecolor = 'lime', alpha = 0.6, linewidth = 2)
                    else:
                        im = ax.bar(x=np.deg2rad(pa), height = height, width=theta_width, bottom=vel, color = color, alpha = 0.7)            
                
            if bellow_cont == 'on':             
                for data_line in reg: #plot data bellow continuum level with red outline
                    if data_line['p'] < cont_p:
                        pa = data_line['pa']; pa_error = data_line['paerr']; vel = data_line['vel']
                        theta_width = np.deg2rad((pa + pa_error) - (pa - pa_error)) #width of line is error of PA
                        height = get_velocity(rest_wave, rest_wave+binsize) #thickness of line is one binsize converted to velocity
                        if data_line['flx'] == fmin: #plot data corresponding to flux min with lime green outline 
                            im = ax.bar(x=np.deg2rad(pa), height = height, width=theta_width, bottom=vel, color = color, edgecolor = 'lime', alpha = 0.6, linewidth = 2)
                        else:
                            im = ax.bar(x=np.deg2rad(pa), height = height, width=theta_width, bottom=vel, color = color, edgecolor = 'red', alpha = 0.6, linewidth = 2)
               
    else:
        #if regions are isolated by hand and specified by [blue_end, red_end] wavelength pairs in isolate use those
        for count, rest_wave in enumerate(lines):
            bluelim_wave = isolate[count][0]
            redlim_wave = isolate[count][1]
            reg = get_lines(epoch_data, bluelim_wave, redlim_wave) #get data for region 
            reg = get_vel_column(reg, rest_wave) #get velocity column 
            reg = get_pa_column(reg, adjust_pa) #get pa and pa_err column ('on' specifies adjust PA = PA +180)
            color = colors[count]
            
            pmax = np.nanmax(reg['p'])
            fmin = np.nanmin(reg['flx'])
            if bellow_cont == 'off':
                for data_line in reg: #plot data above continuum level
                    pa = data_line['pa']; pa_error = data_line['paerr']; vel = data_line['vel']
                    theta_width = np.deg2rad((pa + pa_error) - (pa - pa_error)) #width of line is error of PA
                    if theta_width < .17:
                        theta_width = np.deg2rad(10)
                    height = get_velocity(rest_wave, rest_wave+binsize) #thickness of line is one binsize converted to velocity
                    #if data_line['p'] == pmax: #plot data with pmax black outline 
                        #im = ax.bar(x=np.deg2rad(pa), height = height, width=theta_width, bottom=vel, color = color, edgecolor = 'k')
                        #print(pa, pa_error)
                    #if data_line['flx'] == fmin: #plot data corresponding to flux min with lime green outline and translucent so it does not come out darker 
                        #im = ax.bar(x=np.deg2rad(pa), height = height, width=theta_width, bottom=vel, color = color, edgecolor = 'lime', alpha = 0.6, linewidth = 2)
                    #else:
                    im = ax.bar(x=np.deg2rad(pa), height = height, width=theta_width, bottom=vel, color = color, alpha = .6)#plot other points as solid color            
                        #print(pa, pa_error)
            if bellow_cont == 'on':             
                for data_line in reg: #plot data bellow continuum level as translucent 
                    if data_line['p'] < cont_p:
                        pa = data_line['pa']; pa_error = data_line['paerr']; vel = data_line['vel']
                        theta_width = np.deg2rad((pa + pa_error) - (pa - pa_error)) #width of line is error of PA
                        if theta_width < .17:
                            theta_width = np.deg2rad(10)
                        height = get_velocity(rest_wave, rest_wave+binsize) #thickness of line is one binsize converted to velocity
                        #if data_line['flx'] == fmin: #plot data corresponding to flux min with lime green outline 
                            #im = ax.bar(x=np.deg2rad(pa), height = height, width=theta_width, bottom=vel, color = color, edgecolor = 'lime', alpha = 0.6, linewidth = 2)
                        #else:
                        im = ax.bar(x=np.deg2rad(pa), height = height, width=theta_width, bottom=vel, color = color, alpha = .6, edgecolor = 'r')#plot points bellow cont as translucent 
                        #print(pa, pa_error)
            
            if bellow_cont == 'cs':   #for data after continuum subtraction (cs)        
                for data_line in reg: 
                    pa = data_line['pa']; pa_error = data_line['paerr']; vel = data_line['vel']
                    #print(pa, pa_error)
                    theta_width = np.deg2rad((pa + pa_error) - (pa - pa_error)) #width of line is error of PA
                    if theta_width < .17:
                        theta_width = np.deg2rad(10)
                    height = get_velocity(rest_wave, rest_wave+binsize) #thickness of line is one binsize converted to velocity
                    if data_line['p'] == pmax: #plot data with pmax black outline and solid color
                        #im = ax.bar(x=np.deg2rad(pa), height = height, width=theta_width, bottom=vel, color = color, edgecolor = color, zorder = 2)
                        im = ax.bar(x=np.deg2rad(pa), height = height, width=theta_width, bottom=vel, color = color, edgecolor = 'k', zorder = 2)
                    else:
                        if np.rad2deg(theta_width) < 25:
                            im = ax.bar(x=np.deg2rad(pa), height = height, width=theta_width, bottom=vel, color = color, alpha = .6)#plot all other points as translucent 
                    
    return(im)   

def ppp_val(data, obs_wave, rest_wave, axs, color, binsize = 15):
    """Polar Plot Point Value:
    To plot a single point on polar plot, must use with plot_polar and define axis first.
    obs_wave should be the wavelength associated with a min or max from get_min or get_max.
    Returns polarization value."""
    line = get_lines(data, obs_wave-5, obs_wave+5); line = get_vel_column(line, rest_wave); line = get_pa_column(line, 'on')
    pa = line['pa']; pa_error = line['paerr']; vel = line['vel']
    theta_width = np.deg2rad((pa + pa_error) - (pa - pa_error)) #width of line is error of PA
    height = get_velocity(rest_wave, rest_wave+binsize) #thickness of line is one binsize converted to velocity
    axs.bar(x=np.deg2rad(pa), height = height, width=theta_width, bottom=vel, color = color)
    return(line['p'], line['perr'])      

from skimage.transform import resize
def make_gif(folder, file_list, movie_name, slide_time):
    """makes a gif out of files provided in the file_list. 
    the Gif is saved under the movie_name.GIF (make sure to include .GIF at the end of the movie name) 
    and the slide_time specifies how long to spend on each image"""
    
    images = []

    for filename in file_list:
        images.append(imageio.imread(folder+filename))
    imageio.mimsave(str(movie_name), images, duration = slide_time)
    return(print("GIF saved"))

