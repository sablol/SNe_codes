# -*- coding: utf-8 -*-
"""
Created on Sat Jul 31 09:44:58 2021

@author: sabri
"""
from get_data import *
from Binning_FITS import *
#from plot_data import *
#from astropy.io import ascii 
from plot_codes import *
from astropy.table import Table
import matplotlib.pyplot as plt
from find_lines import *
import numpy as np

#open, combine and bin data
unbin1 = get_txtFITS("Epoch_1/all_comb", "sn2012au_allcomb.qsum.txt", "sn2012au_allcomb.q.txt", "sn2012au_allcomb.qsig.txt", "sn2012au_allcomb.flx.txt", "sn2012au_allcomb.u.txt", "sn2012au_allcomb.usig.txt")
epoch1_20 = Bin_data(unbin1, 4, 20)
data2 = get_txtFITS("Epoch_2", "sn2012au_allcomb.qsum.txt", "sn2012au_allcomb.q.txt", "sn2012au_allcomb.qsig.txt", "sn2012au_allcomb.flx.txt", "sn2012au_allcomb.u.txt", "sn2012au_allcomb.usig.txt")
epoch2_20 = Bin_data(data2, 2.5, 20)
data3 = get_txtFITS("Epoch_3", "sn2012au_allcomb.qsum.txt", "sn2012au_allcomb.q.txt", "sn2012au_allcomb.qsig.txt", "sn2012au_allcomb.flx.txt", "sn2012au_allcomb.u.txt", "sn2012au_allcomb.usig.txt")
epoch3_20 = Bin_data(data3, 4, 20)
data4 = get_txtFITS("Epoch_4", "sn2012au_allcomb.qsum.txt", "sn2012au_allcomb.q.txt", "sn2012au_allcomb.qsig.txt", "sn2012au_allcomb.flx.FIXED.txt", "sn2012au_allcomb.u.txt", "sn2012au_allcomb.usig.txt")
epoch4_20 = Bin_data(data4, 4, 20); epoch4_20 = epoch4_20[10:]
data5 = get_txtFITS("Epoch_5", "sn2012au_allcomb.qsum.txt", "sn2012au_allcomb.q.txt", "sn2012au_allcomb.qsig.txt", "sn2012au_allcomb.flx.txt", "sn2012au_allcomb.u.txt", "sn2012au_allcomb.usig.txt")
epoch5_20 = Bin_data(data5, 4, 20)
data6 = get_fits("Epoch_6")
data6.rename_column('usum', 'flx')
epoch6_20 = Bin_data(data6, 2.5, 20)

#-------ALL DATA----------------
"""
#comparing %p calcuations------------------------------------------------------------
oldp = P_tot(epoch1_20['q'], epoch1_20['u'])
newp = P_debi(epoch1_20['q'], epoch1_20['u'], epoch1_20['qerr'], epoch1_20['uerr'])
plt.step(epoch1_20['wave'], oldp, label = ('OG'))
plt.step(epoch1_20['wave'], newp, label = ('debias'))
plt.xlabel("wavelength"); plt.ylabel("total %P")
plt.legend()
"""
"""
#pannel plots-------------------------------------------------------------------------------------
plot_pannels(unbin1, epoch1_20, 10**-4.5, 10**-4.5, title = "Days 0-7 Binned to 20 Angstroms")
plot_pannels(data2, epoch2_20, 10**-4.6, 10**-4.7, title = "Day 26 Binned to 20 Angstroms")
plot_pannels(data3, epoch3_20, 10**-4, 10**-4.4, title = "Days 35-40 Binned to 20 Angstroms")
plot_pannels(data4, epoch4_20, 10**-3.7, 10**-3.7, title = "Days 57-67 Binned to 20 Angstroms")
plot_pannels(data5, epoch5_20, 10**-3.2, 10**-3.2, title = "Days 85-90 Binned to 20 Angstroms")
plot_pannels(data6, epoch6_20, 10**-3.8, 10**-3.8, 'off', title = "Day 295 Binned to 20 Angstroms")
"""
"""
#QU whole epoch plots-------------------------------------------------------------------------------
plot_QU(epoch1_20, "Days0-7_binned20")
plot_QU(epoch2_20, "Epoch 2: Day 26")
#plot_QU(epoch3_20, "Days35-40_binned20")
#plot_QU(epoch4_20, "Days57-67_binned20")
#plot_QU(epoch5_20, "Days85-90_binned20")
#plot_QU(epoch6_20, "Day295_binned20")
"""
"""
#all epoch total P and flux and lines plots---------------------------------------------------------------
all_epochs20 = [epoch1_20, epoch2_20, epoch3_20, epoch4_20, epoch5_20, epoch6_20]
all_epochs_unbin = [unbin1, data2, data3, data4, data5, data6]
early = [epoch1_20, epoch2_20, epoch3_20]; early_unbin = [unbin1, data2, data3]
late = [epoch4_20, epoch5_20, epoch6_20]

flx_adjust = [10**-5, 2*10**-6, 3*10**-5, 5*10**-5, 7*10**-5, 1.6*10**-5 ]
flx_adjust_early = [.00005, .00002, .0002]
flx_adjust_late = [.0003, .0004, .00015]
polflx_adjust = [10**-5, 2*10**-6, 4*10**-5, 6*10**-5, 8*10**-5, 3*10**-5]

epoch_names = ['days 0-7', 'day 26', 'days 35-40', 'days 57-67', 'days 85-90', 'day 295']
early_epochs = ['days 0-7', 'day 26', 'days 35-40']
late_epochs = ['days 57-67', 'days 85-90', 'day 295']
"""
"""
#all data plots (messy but good initial look tool)
#plot_all_data(all_epochs_unbin, all_epochs20, flx_adjust, polflx_adjust, epoch_names)

#all flux data stacked
flx_adjust = [10**-4, 3*10**-5, 3*10**-4, 5*10**-4, 7*10**-4, 1.6*10**-4]
epoch_names = ['days 0-7', 'day 26', 'days 35-40', 'days 57-67', 'days 85-90', 'day 295']
plot_all_flx(all_epochs_unbin, flx_adjust, epoch_names, 6)
"""
"""
##### Line IDs
He_lines = [4471,4922, 5016, 5876, 6678, 7065]
He_labels = ["4471","4922", "5016", "5876", "6678", "7065"]
He_colors = ['blue', "blue", "blue", 'blue', 'blue', 'blue']
He_regions = [[5400,5800], [6300, 6400], [6800, 6900]]
lines = [4471, 4924, 5876, 6678, 7065, 4571, 5169, 6300, 6364, 4959, 5007, 6355, 6800, 7600]
line_c = ['b', 'g', 'b', 'b', 'b', 'm', 'g', 'y', 'y', 'y', 'y', 'r', 'k', 'k']
line_lab = ['He I', 'Fe II', 'He+NaID', '', '', 'Mg', 'Fe II', 'OI', '', '', 'OIII', 'SiII', '$\\bigoplus$', '$\\bigoplus$']
penday_lines = [4571, 5000, 5363, 5876, 5535, 5890, 6355, 6678, 7065, 7774]
penday_labels = ['Mg I 4571', 'Fe Triplet', 'Fe II 5363', '', 'Fe II 5535', 'Na ID 5890', 'Ha +SiII 6355', 'He I 6678', 'He I 7065', 'O I 7774']
penday_colors = ['m','g', 'g', 'b', 'g', 'r', 'c', 'b', 'b', 'y']
p_early_lines = [4471, 5876, 6678, 7065, 6355, 4549, 4925, 5018, 5169]
p_early_labels = ['He 4471', 'He 5876', 'He 6678', 'He 7065', 'Si 6355', 'Fe II 4549', 'Fe II Triplet', '', 'Fe 5169']
p_early_colors = ['b', 'b', 'b', 'b', 'c', 'g', 'g', 'g', 'g']
p_early_regions = [[4300, 4500], [4800, 5000], [5600, 5800], [6000, 6100], [6300, 6400], [6600, 6700]]
p_late_lines = [4571, 5000, 5363, 5535, 5890, 6350]
p_late_labels = ['Mg 4571', 'Fe Triplet', 'Fe 5336', 'Fe 5535', 'Na ID 5890', 'O I doublet']
p_late_colors = ['m', 'g', 'g', 'g', 'r', 'y']

plot_all_epochs(all_epochs20, flx_adjust, epoch_names, xmin = 3700, title = "", line_list = He_lines, color_list = He_colors , label_list = He_labels) 
#plot_all_epochs(all_epochs20, flx_adjust, epoch_names, xmin = 3700, line_list = penday_lines, label_list = penday_labels, color_list = penday_colors, title = "") 
#im = early_spec = plot_all_epochs(early, flx_adjust_early, early_epochs, xmin = 3700, line_list = He_lines, label_list= He_labels, color_list= He_colors, title = " ")
#im.axvline(x = 5582, ymin = .85, ymax = .98, c = 'k'); im.axvline(x = 6345, ymin = .85, ymax = .90, c = 'k');im.axvline(x = 6850, ymin = .82, ymax = .87, c = 'k')
#im.axvline(x = 5650, ymin = .44, ymax = .61, c = 'k'); im.axvline(x = 6440, ymin = .52, ymax = .57, c = 'k'); im.axvline(x = 6850, ymin = .47, ymax = .53, c = 'k')
#im.axvline(x = 5680, ymin = .15, ymax = .28, c = 'k'); im.axvline(x = 6530, ymin = .21, ymax = .26, c = 'k'); im.axvline(x = 6850, ymin = .12, ymax = .20, c = 'k')
#late_spec = plot_all_epochs(late, flx_adjust_late, late_epochs, xmin = 3700, line_list = p_late_lines, label_list= p_late_labels, color_list= p_late_colors, title = " ")
"""
"""
#------------CONTINUUM REGIONS-------------------
cont_flx_adjust = [10**-5, 3*10**-6, 3*10**-5, 5*10**-5, 7*10**-5, 1.6*10**-5]
continuum_regions = [[5200, 5400], [7000, 7300]]
plot_all_epochs(all_epochs20, cont_flx_adjust, epoch_names, xmin = 3700, regions_list = continuum_regions, title = "") 
"""
"""
#continuum region averages Q/U plot-------------------------------------------------------------------------------
all_epochs = [epoch1_20, epoch2_20, epoch3_20, epoch4_20, epoch5_20, epoch6_20]
cont1 = get_aves(all_epochs, [5200, 5400])
epochs1_5 = [epoch1_20, epoch2_20, epoch3_20, epoch4_20, epoch5_20]
cont2 = get_aves(all_epochs, [7000, 7300])
cont1_PAs = PA(cont1['q'], cont1['u']); cont1_PAerr = PA_err(cont1['p'], cont1['q'], cont1['u'], cont1['qerr'], cont1['uerr']) 
cont2_PAs = PA(cont2['q'], cont2['u']); cont2_PAerr = PA_err(cont2['p'], cont2['q'], cont2['u'], cont2['qerr'], cont2['uerr'])   

#our ------ISP estimate----- from weighted average of epoch 5 continuum regions
cont1[5] = cont2[4] #add second region epoch 5 data to table
cont1['wave'][0] = 0 #reset all other epoch wavelengths so they are not considered in average 
cont1['wave'][1] = 0
cont1['wave'][2] = 0
cont1['wave'][3] = 0
isp = get_aves([cont1], [5293, 7144]) #take average of two regions
isp_q = isp['q']; isp_u = isp['u']; #print(isp_q, isp_u)
print(isp)
print(PA(isp_q, isp_u))
print(PA_err(isp['p'], isp_q, isp_u, isp['qerr'], isp['uerr']))


#cont1 sets up axis
plot_QU(cont1, 'SN 2012au Continuum Regions 5100-5500 & 7000-7300', size = 0.02, epoch_labels = ['1', '2', '3', '4', '5', '6'])
plt.plot(cont1['q'], cont1['u'], c = 'b', linewidth = 1, label= '5100-5500Å')
#overplot cont2
plt.plot(cont2['q'], cont2['u'], c = 'r', linewidth = 1, label= '7000-7300Å')
plt.scatter(cont2['q'], cont2['u'], marker = 'o', s=20)
plt.errorbar(cont2['q'], cont2['u'], yerr=cont2['uerr'], xerr=cont2['qerr'], ecolor="grey", hold=True, fmt="none")
epoch1_5labels = ['1', '2', '3', '4', '5']
for i, txt in enumerate(epoch1_5labels):
    plt.annotate(txt, xy = (cont2['q'][i]+0.0005, cont2['u'][i]), fontsize = 14) #label each data point according to epoch_labels provided     
#best fit line from epoch 1
Q1 = [x for x in epoch1_20['q'] if np.isnan(x) == False] #call data and ignore nan values
U1 = [x for x in epoch1_20['u'] if np.isnan(x) == False]
Qsig1 = [x for x in epoch1_20['qerr'] if np.isnan(x) == False]
Usig1 = [x for x in epoch1_20['uerr'] if np.isnan(x) == False]    
bfit = best_fit(Q1, U1, Qsig1, Usig1)
plt.plot(bfit[0], bfit[1], 'c', label = 'best fit epoch1')
plt.legend()


#continuum line fit estimate from 5876 all epoch Q/U plots---------------------------------
fig, ax = plt.subplots(2, 3, figsize = (30, 16)); fig.subplots_adjust(right=0.8)
fig, ax = plt.subplots(2, 3, figsize = (30, 16), sharex = True, sharey = True); fig.subplots_adjust(right = .85,  wspace=-.33, hspace=0.05)
fig.text(0.5, 0.08, '% Q Polarization', ha='center'); fig.text(0.15, 0.5, '% U Polarization', va='center', rotation='vertical') #set shared axis labels 

#local continuum region data
all_epochs = [epoch1_20, epoch2_20, epoch3_20, epoch4_20, epoch5_20, epoch6_20[0:-12]]
cont1 = get_aves(all_epochs, [5100, 5500]) 
qcont = cont1['q']; ucont = cont1['u']; qcerr = cont1['qerr']; ucerr = cont1['uerr']
ax[0,0].plot(qcont[0], ucont[0], marker = '*', c = 'k', markersize = 20, label = 'continuum'); ax[0,1].plot(qcont[1], ucont[1], marker = '*', c = 'k', markersize = 20, label = 'continuum')
ax[0,2].plot(qcont[2], ucont[2], marker = '*', c = 'k', markersize = 20, label = 'continuum'); ax[1,0].plot(qcont[3], ucont[3], marker = '*', c = 'k', markersize = 20)
ax[1,1].plot(qcont[4], ucont[4], marker = '*', c = 'k', markersize = 20); ax[1,2].plot(qcont[5], ucont[5], marker = '*', c = 'k', markersize = 20)
#get data for area best fit; take out extruding points for best fit; get data for whole are points
e1bfit = get_lines(epoch1_20, 5530, 6210); e1bfit.remove_rows(slice(3, 12)); e1 = get_lines(epoch1_20, 5530, 6210)
e2bfit = get_lines(epoch2_20, 5530, 6210); e2bfit.remove_rows(slice(3, 12)); e2 = get_lines(epoch2_20, 5530, 6210)
e3bfit = get_lines(epoch3_20, 5530, 6210); e3bfit.remove_rows(slice(3, 12)); e3 = get_lines(epoch3_20, 5530, 6210)
e4bfit = get_lines(epoch4_20, 5530, 6210); e4bfit.remove_rows(slice(3, 12)); e4 = get_lines(epoch4_20, 5530, 6210)
e5bfit = get_lines(epoch5_20, 5530, 6210); e5bfit.remove_rows(slice(3, 12)); e5 = get_lines(epoch5_20, 5530, 6210)
e6bfit = get_lines(epoch6_20, 5530, 6210); e6bfit.remove_rows(slice(3, 12)); e6 = get_lines(epoch6_20, 5530, 6210)
extruding_points = e1['wave'][3:12]
#calculate best fit for data without extruding points 
bfit1 =  
bfit2 = best_fit(e2bfit['q'], e2bfit['u'], e2bfit['qerr'], e2bfit['uerr'])
bfit3 = best_fit(e3bfit['q'], e3bfit['u'], e3bfit['qerr'], e3bfit['uerr'])
bfit4 = best_fit(e4bfit['q'], e4bfit['u'], e4bfit['qerr'], e4bfit['uerr'])
bfit5 = best_fit(e5bfit['q'], e5bfit['u'], e5bfit['qebest_fit(e1bfit['q'], e1bfit['u'], e1bfit['qerr'], e1bfit['uerr'])rr'], e5bfit['uerr'])
bfit6 = best_fit(e6bfit['q'], e6bfit['u'], e6bfit['qerr'], e6bfit['uerr'])
#Plot QU for all points in area and best fit without extruding points 
epoch1 = QU(e1, ax[0,0], bfit = bfit1); epoch2 = QU(e2, ax[0,1], bfit = bfit2); epoch3 = QU(e3, ax[0,2], bfit = bfit3)
epoch4 = QU(e4, ax[1,0], bfit = bfit4); epoch5 = QU(e5, ax[1,1], bfit = bfit5); epoch6 = QU(e6, ax[1,2], bfit = bfit6)
cbar_ax = fig.add_axes([0.81, 0.18, 0.04, 0.68])
cbar = fig.colorbar(epoch1, cax=cbar_ax).set_label('(Å)', rotation= 0, fontsize = 20)
ax[0,0].text(1.6, 3.4, 'Days 0-7'); ax[0,1].text(1.6, 3.4, 'Day 26');ax[0,2].text(1.6, 3.4, 'Days 35-40')
ax[1,0].text(1.6, 3.4, 'Days 57-67'); ax[1,1].text(1.6, 3.4, 'Days 85-90');ax[1,2].text(1.6, 3.4, 'Day 295')

#plot rest wavelength with special marker:
ax[0,0].plot(e1['q'][16], e1['u'][16], marker = 'X', c = 'lime', markersize = 20, label = "rest\nwavelength")
ax[0,1].plot(e2['q'][16], e2['u'][16], marker = 'X', c = 'lime', markersize = 20, label = "rest λ")
ax[0,2].plot(e3['q'][17], e3['u'][17], marker = 'X', c = 'lime', markersize = 20, label = "rest λ")
ax[1,0].plot(e4['q'][16], e4['u'][16], marker = 'X', c = 'lime', markersize = 20)
ax[1,1].plot(e5['q'][17], e5['u'][17], marker = 'X', c = 'lime', markersize = 20)
ax[1,2].plot(e6['q'][17], e6['u'][17], marker = 'X', c = 'lime', markersize = 20)
"""
"""
#all epochs - points extruding from 5876 line region best fit highlighted------------------------------------------------------------------------------------------
all_epochs20 = [epoch1_20, epoch2_20, epoch3_20, epoch4_20, epoch5_20[5:], epoch6_20]
flx_adjust = [10**-5, 10**-5.7, 10**-4.5, 10**-4.3, 10**-4, 10**-4.7]#[.00005, .00002, .0002, .0003, .0004, .00015]
epoch_names = ['days 0-7', 'day 26', 'days 35-40', 'days 57-67', 'days 85-90', 'day 295']
plot_all_epochs(all_epochs20, flx_adjust, epoch_names, xmin = 3700, regions_list = [[5200, 5400], [7000, 7300]], title = "SN 2012au")
"""
#--------------------ISP from all estimates-------------
"""
#####Get late time continuum region averages
#from astropy.table import Table, hstack, vstack
#cont1 = get_aves([epoch5_20, epoch6_20], [5200, 5400])
#cont2 = get_aves([epoch5_20], [7000, 7300])
#e5e6cont = vstack([cont1, cont2]);
#cont_ave = get_aves([e5e6cont], [5291, 7145])
#print(cont_ave['q'], cont_ave['u'], cont_ave['qerr'], cont_ave['uerr'])
#####set up figure
fig, ax = plt.subplots(1, figsize = (10, 10), sharex = True, sharey = True); fig.subplots_adjust(right = .85,  wspace=-.33, hspace=0.05)
fig.text(0.5, 0.08, '% Q Polarization', ha='center'); fig.text(0.05, 0.5, '% U Polarization', va='center', rotation='vertical') #set shared axis labels 
#####ISP sing point esitmates [Pandey probe star fit, my probe star, my cont. epoch 5, cont epoch 5 + 6]
qisp = [-0.1407, 0.1099, 0.17, .22]; uisp = [0.1819, 0.2433, 0.18, .19]; qerrisp = [0.01, 0.036, 0.037, .017]; uerrisp = [0.01, 0.036, 0.042, .021]
ax.errorbar(qisp[0], uisp[0], xerr = qerrisp[0], yerr = uerrisp[0], c = 'b', label = 'Pandey et al probe stars')
ax.errorbar(qisp[1], uisp[1], xerr = qerrisp[1], yerr = uerrisp[1], c = 'g', label = 'HD 112142')
#ax.errorbar(qisp[2], uisp[2], xerr = qerrisp[2], yerr = uerrisp[2], c = 'purple', label = 'Epoch 5 continuum ave.')
#ax.errorbar(qisp[3], uisp[3], xerr = qerrisp[3], yerr = uerrisp[3], c = 'hotpink', label = 'all Epoch 5+ blue epoch 6 continuum ave.')

#####ISP 9xE(B-V) circle 
ax.add_patch(plt.Circle((0, 0), radius = .32, edgecolor='red', facecolor = 'none', label = 'Serkwoski-Reddening'))
"""
"""
#continuum region by epoch continuum_regions = [[5200, 5400], [7000, 7300]]
e1bfit1 = get_lines(epoch1_20, 5200, 5400)
bfit = best_fit(e1bfit1['q'], e1bfit1['u'], e1bfit1['qerr'], e1bfit1['uerr'])
ax.plot(bfit[0], bfit[1], 'b', label = 'best fit epoch1 (5200-5400)')
e1bfit2 = get_lines(epoch1_20, 7000, 7300)
bfit = best_fit(e1bfit2['q'], e1bfit2['u'], e1bfit2['qerr'], e1bfit2['uerr'])
ax.plot(bfit[0], bfit[1], 'b:', label = 'best fit epoch1 (7000-7300)')

e2bfit1 = get_lines(epoch2_20, 5200, 5400)
bfit = best_fit(e2bfit1['q'], e2bfit1['u'], e2bfit1['qerr'], e2bfit1['uerr'])
ax.plot(bfit[0], bfit[1], 'y', label = 'best fit epoch2 (5200-5400)')
e2bfit2 = get_lines(epoch2_20, 7000, 7300)
bfit = best_fit(e2bfit2['q'], e2bfit2['u'], e2bfit2['qerr'], e2bfit2['uerr'])
ax.plot(bfit[0], bfit[1], 'y:', label = 'best fit epoch2 (7000-7300)')

e3bfit1 = get_lines(epoch3_20, 5200, 5400)
bfit = best_fit(e3bfit1['q'], e3bfit1['u'], e3bfit1['qerr'], e3bfit1['uerr'])
ax.plot(bfit[0], bfit[1], 'g', label = 'best fit epoch3 (5200-5400)')
e3bfit2 = get_lines(epoch3_20, 7000, 7300)
bfit = best_fit(e3bfit2['q'], e3bfit2['u'], e3bfit2['qerr'], e3bfit2['uerr'])
ax.plot(bfit[0], bfit[1], 'g:', label = 'best fit epoch3 (7000-7300)')

e4bfit1 = get_lines(epoch4_20, 5200, 5400)
bfit = best_fit(e4bfit1['q'], e4bfit1['u'], e4bfit1['qerr'], e4bfit1['uerr'])
ax.plot(bfit[0], bfit[1], 'r', label = 'best fit epoch4 (5200-5400)')
e4bfit2 = get_lines(epoch4_20, 7000, 7300)
bfit = best_fit(e4bfit2['q'], e4bfit2['u'], e4bfit2['qerr'], e4bfit2['uerr'])
ax.plot(bfit[0], bfit[1], 'r:', label = 'best fit epoch4 (7000-7300)')

e5bfit1 = get_lines(epoch5_20, 5200, 5400)
bfit = best_fit(e5bfit1['q'], e5bfit1['u'], e5bfit1['qerr'], e5bfit1['uerr'])
ax.plot(bfit[0], bfit[1], 'purple', label = 'best fit epoch5 (5200-5400)')
e5bfit2 = get_lines(epoch5_20, 7000, 7300)
bfit = best_fit(e5bfit2['q'], e5bfit2['u'], e5bfit2['qerr'], e5bfit2['uerr'])
ax.plot(bfit[0], bfit[1], 'purple', linestyle = ':', label = 'best fit epoch5 (7000-7300)')

e6bfit1 = get_lines(epoch6_20, 5200, 5400)
bfit = best_fit(e6bfit1['q'], e6bfit1['u'], e6bfit1['qerr'], e6bfit1['uerr'])
ax.plot(bfit[0], bfit[1], 'brown', label = 'best fit epoch6 (5200-5400)')
e6bfit2 = get_lines(epoch6_20, 7000, 7300)
bfit = best_fit(e6bfit2['q'], e6bfit2['u'], e6bfit2['qerr'], e6bfit2['uerr'])
ax.plot(bfit[0], bfit[1], 'brown', linestyle = ':', label = 'best fit epoch6 (7000-7300)')
"""
"""
e1bfit = get_lines(epoch1_20, 5530, 6210); e1bfit.remove_rows(slice(3, 12)); 
e2bfit = get_lines(epoch2_20, 5530, 6210); e2bfit.remove_rows(slice(3, 12));
e1e2bfit = hstack([e1bfit, e2bfit])
bfit = best_fit(e1e2bfit['q'], e1e2bfit['u'], e1e2bfit['qerr'], e1e2bfit['uerr'])
ax.plot(bfit[0], bfit[1], 'c', label = 'best fit epochs1+2 (no 5876 region)')

e1bfit = get_lines(epoch1_20, 5530, 6210)#; e1bfit.remove_rows(slice(3, 12)); 
e2bfit = get_lines(epoch2_20, 5530, 6210)#; e2bfit.remove_rows(slice(3, 12));
e1e2bfit = hstack([e1bfit, e2bfit])
bfit = best_fit(e1e2bfit['q'], e1e2bfit['u'], e1e2bfit['qerr'], e1e2bfit['uerr'])
ax.plot(bfit[0], bfit[1], 'm', label = 'best fit all epochs1+2')

e1bfit = get_lines(epoch1_20, 5530, 6210)#; e1bfit.remove_rows(slice(3, 12)); 
bfit = best_fit(e1bfit['q'], e1bfit['u'], e1bfit['qerr'], e1bfit['uerr'])
ax.plot(bfit[0], bfit[1], 'y', label = 'best fit all epoch1')

e2bfit = get_lines(epoch2_20, 5530, 6210); e2bfit.remove_rows(slice(3, 12)); 
bfit = best_fit(e2bfit['q'], e2bfit['u'], e2bfit['qerr'], e2bfit['uerr'])
ax.plot(bfit[0], bfit[1], 'darkorange', label = 'best fit epoch2 (no 5876region)')
"""
"""
#####continuum best fit lines
all_epochs = [epoch1_20, epoch2_20, epoch3_20, epoch4_20, epoch5_20, epoch6_20]
cont1 = get_aves(all_epochs, [5200, 5400])
cont2 = get_aves(all_epochs, [7000, 7300])
cont_tot = hstack([cont1, cont2])   
bfit = best_fit(cont_tot['q'], cont_tot['u'], cont_tot['qerr'], cont_tot['uerr'])
ax.plot(bfit[0], bfit[1], 'k', label = 'allepochs both regions continuum best fit');
"""
"""
cont1_epochs = [epoch3_20, epoch4_20, epoch5_20, epoch6_20]
cont1 = get_aves(cont1_epochs, [5200, 5400])
bfit = best_fit(cont1['q'], cont1['u'], cont1['qerr'], cont1['uerr'])
ax.plot(bfit[0], bfit[1], 'c', label = 'epochs3,4,5,6 5200-5400 continuum best fit');
"""
"""
cont1_epochs = [epoch3_20, epoch4_20, epoch5_20, epoch6_20]
cont1 = get_aves(cont1_epochs, [5100, 5500])
bfit = best_fit(cont1['q'], cont1['u'], cont1['qerr'], cont1['uerr'])
ax.plot(bfit[0], bfit[1], 'darkorange', label = 'epochs3,4,5,6 5100-5500 continuum best fit');
Reilly16ISPB = get_aves([cont1], [5200, 5300])
ax.errorbar(Reilly16ISPB['q'], Reilly16ISPB['u'], xerr = Reilly16ISPB['qerr'], yerr = Reilly16ISPB['uerr'], c = 'c', label = 'Reilly et al ISP (5100 - 5500) region ave')
"""
"""
cont2_epochs = [epoch3_20, epoch4_20, epoch5_20]
cont2 = get_aves(cont2_epochs, [7000, 7300])
bfit = best_fit(cont2['q'], cont2['u'], cont2['qerr'], cont2['uerr'])
ax.plot(bfit[0], bfit[1], 'g', label = 'epochs3,4,5 7000-7300 continuum best fit');
#####epoch 2 He 5876 data for comparison 
#e2 = plot_line_velspace(epoch2_20, 5876, 20000, 10**-5.7, "Epoch2(Day26)")
#epoch2 = QU(e2, ax)
"""
"""
#####Reilly 2016 cont. regions best fit lines
all_epochs = [epoch1_20, epoch3_20, epoch4_20, epoch5_20]
cont1 = get_aves(all_epochs, [6050, 6225])
cont2 = get_aves(all_epochs, [7100, 7500])
cont_tot = hstack([cont1, cont2])   
bfit = best_fit(cont_tot['q'], cont_tot['u'], cont_tot['qerr'], cont_tot['uerr'])
ax.plot(bfit[0], bfit[1], 'y', label = 'epochs 1,3,4,5 Reilly cont regions 5050-5225, 7100-7500')
"""
"""
all_epochs = [epoch1_20, epoch3_20, epoch4_20, epoch5_20]
cont2 = get_aves(all_epochs, [7100, 7500])   
bfit = best_fit(cont2['q'], cont2['u'], cont2['qerr'], cont2['uerr'])
ax.plot(bfit[0], bfit[1], 'm', label = 'epochs 1,3,4,5 Reilly cont region 7100-7500');
"""
"""
#####set axis 
ax.set_title('Continuum region best fit by epoch with ISP estimates')
ax.axhline(0, color = 'k', linewidth = 1)
ax.axvline(0, color = 'k', linewidth=1)
ax.legend()
ax.axis('square')
ax.set(xlim=(-1.5, 1.5), ylim=(-1.5, 1.5))
"""
#-----------HE 5876 LINE------------------
"""
#all epochs 5876 line region highlighted------------------------------------------------------------------------------------------
all_epochs20 = [epoch1_20, epoch2_20, epoch3_20, epoch4_20, epoch5_20, epoch6_20[0:-12]]
flx_adjust = [10**-5, 10**-5.7, 10**-4.5, 10**-4.3, 10**-4, 10**-4.7]
epoch_names = ['days 0-7', 'day 26', 'days 35-40', 'days 57-67', 'days 85-90', 'day 295']
#print(get_obswave(5876, 20000))
plot_all_epochs(all_epochs20, flx_adjust, epoch_names, xmin = 5500, xmax = 6300, title = "SN 2012au He I +/-20000km/s", line_list = [5876], color_list = ['b'], label_list = ['He I 5876'])
"""
"""
all_epochs20 = [epoch1_20, epoch2_20, epoch3_20, epoch4_20, epoch5_20, epoch6_20[0:-12]]
flx_adjust = [10**-5, 10**-5.6, 10**-4.3, 10**-4.1, 10**-4.1, 10**-4.5]
epoch_names = ['days 0-7', 'day 26', 'days 35-40', 'days 57-67', 'days 85-90', 'day 295']
print(get_obswave(5876, 17000))
plot_all_epochs(all_epochs20, flx_adjust, epoch_names, xmin = 5540, xmax = 6200, line_list = [5876], color_list = ['k'], label_list = ['He I 5876'], polflux= 'banana')
"""
"""
#Velocity Space plots - line 5876 -------------------------------------------------------------------
e1 = plot_line_velspace(epoch1_20, 5876, 15000, 10**-5,"Epoch1(Day0-7)")
e2 = plot_line_velspace(epoch2_20, 5876, 15000, 10**-5.7, "Epoch2(Day26)")
e3 = plot_line_velspace(epoch3_20, 5876, 15000, 10**-4.5, "Epoch3(Day35-40)")
e4 = plot_line_velspace(epoch4_20, 5876, 15000, 10**-4.3, "Epoch4(Day57-67)")
e5 = plot_line_velspace(epoch5_20, 5876, 15000, 10**-4, "Epoch5(Day85-90)")
e6 = plot_line_velspace(epoch6_20, 5876, 15000, 10**-4.7, "Epoch6(Day295)")
#PAer = PA_err(e3['p'], e3['q'], e3['u'], e3['qerr'], e3['uerr'])
#print(np.mean(e3['perr']))
"""
"""
#plot of epochs 1-3 with He 5876 velocity shift labeled ------------------------------------------------------------
epochs1_3_20 = [epoch1_20, epoch2_20, epoch3_20]
flx_adjust = [.00005, .00002, .0002]
epoch_names = ['days 0-7', 'day 26', 'days 35-40']
im = plot_all_epochs(epochs1_3_20, flx_adjust, epoch_names, xmin = 3700, line_list=[4471, 4924, 5876, 6678, 7065], color_list= ['b', 'b', 'b', 'b', 'b'], label_list= ["4471", "4924", "5876", "6678", "7065"]) 
#im.axvline(x = 5582, ymin = .70, ymax = .98, c = 'c'); im.text(5250, 16.5, '-15,000 km/s') #velocities based on flux minimums
#im.axvline(x = 5650, ymin = .35, ymax = .65, c = 'c'); im.text(5250, 11, '-11,530 km/s')
#im.axvline(x = 5665, ymin = .05, ymax = .30, c = 'c'); im.text(5250, 5.5, '-10,765 km/s')
"""
"""
#Single epoch highlight line 5876-------------------------------------------------------------------------------------
He = [5876]; He_lab = ["He I (5876)"]; He_c = ["b"]
#epoch1
plot_p_single(He, He_lab, He_c, unbin1, epoch1_20, 20, 10**-4.5, "Epoch_1")
plt.plot(e1['wave'], e1['p'], c= 'fuchsia', label = "+/- 15000km/s");plt.legend(loc="upper right", fontsize = 12)
plt.savefig("He_5876/Epoch1_binned20_15000kms")
#epoch2
plot_p_single(He, He_lab, He_c, data2, epoch2_20, 20, 10**-4.7, "Epoch_2")
plt.plot(e2['wave'], e2['p'], c= 'fuchsia', label = "+/- 15000km/s");plt.legend(loc="upper right", fontsize = 12)
plt.savefig("He_5876/Epoch2_binned20_15000kms")
#epoch3
plot_p_single(He, He_lab, He_c, data3, epoch3_20, 20, 10**-3.6, "Epoch_3")
plt.plot(e3['wave'], e3['p'], c= 'fuchsia', label = "+/- 15000km/s");plt.legend(loc="upper right", fontsize = 12)
plt.savefig("He_5876/Epoch3_binned20_15000kms")
#epoch4
plot_p_single(He, He_lab, He_c, data4, epoch4_20, 20, 10**-3.5, "Epoch_4")
plt.plot(e4['wave'], e4['p'],c= 'fuchsia', label = "+/- 15000km/s");plt.legend(loc="upper right", fontsize = 12)
plt.savefig("He_5876/Epoch4_binned20_15000kms")
#epoch5
plot_p_single(He, He_lab, He_c, data5, epoch5_20, 20, 10**-3.2, "Epoch_5")
plt.plot(e5['wave'], e5['p'],c= 'fuchsia', label = "+/- 15000km/s");plt.legend(loc="upper right", fontsize = 12)
plt.savefig("He_5876/Epoch5_binned20_15000kms")
#epoch6
plot_p_single(He, He_lab, He_c, data6, epoch6_20, 20, 10**-3.8, "Epoch_6")
plt.plot(e6['wave'], e6['p'], c= 'fuchsia', label = "+/- 15000km/s"); plt.legend(loc="upper right", fontsize = 12)
plt.savefig("He_5876/Epoch6_binned20_15000kms")
"""
"""
#QU all epoch subplots combined - line 5876--------------------------------------------------------------------------------------
import matplotlib.pyplot as plt
fig, ax = plt.subplots(2, 3, figsize = (30, 16), sharex = True, sharey = True); fig.subplots_adjust(right = .85,  wspace=-.33, hspace=0.05)
fig.text(0.5, 0.08, '% Q Polarization', ha='center'); fig.text(0.15, 0.5, '% U Polarization', va='center', rotation='vertical') #set shared axis labels 

#local continuum region data
all_epochs = [epoch1_20, epoch2_20, epoch3_20, epoch4_20, epoch5_20, epoch6_20[0:-12]]
cont1 = get_aves(all_epochs, [5100, 5500]) 
qcont = cont1['q']; ucont = cont1['u']; qcerr = cont1['qerr']; ucerr = cont1['uerr']
cont2 = get_aves(all_epochs, [7000, 7300])
qcont2 = cont2['q']; ucont2 = cont2['u']; qcerr2 = cont2['qerr']; ucerr2 = cont2['uerr']
#local continuum estimates as black star (5100-5500)
ax[0,0].plot(qcont[0], ucont[0], marker = '*', c = 'k', markersize = 20, label = "5100-5500 ave. continuum"); ax[0,1].plot(qcont[1], ucont[1], marker = '*', c = 'k', markersize = 20)
ax[0,2].plot(qcont[2], ucont[2], marker = '*', c = 'k', markersize = 20); ax[1,0].plot(qcont[3], ucont[3], marker = '*', c = 'k', markersize = 20)
ax[1,1].plot(qcont[4], ucont[4], marker = '*', c = 'k', markersize = 20); ax[1,2].plot(qcont[5], ucont[5], marker = '*', c = 'k', markersize = 20)
#red continuum estimate as red star (7000-7300)
ax[0,0].plot(qcont2[0], ucont2[0], marker = '*', c = 'r', markersize = 20, label = "7000-7300 ave. continuum"); ax[0,1].plot(qcont2[1], ucont2[1], marker = '*', c = 'r', markersize = 20)
ax[0,2].plot(qcont2[2], ucont2[2], marker = '*', c = 'r', markersize = 20); ax[1,0].plot(qcont2[3], ucont2[3], marker = '*', c = 'r', markersize = 20)
ax[1,1].plot(qcont2[4], ucont2[4], marker = '*', c = 'r', markersize = 20); ax[1,2].plot(qcont2[5], ucont2[5], marker = '*', c = 'r', markersize = 20)
#Pandey et al. 21 ISP estimate as grey cross
#ax[0,0].plot(0.10133, 0.206476, marker = 'X', c = 'grey', markersize = 20); ax[0,1].plot(0.10133, 0.206476, marker = 'X', c = 'grey', markersize = 20)
#ax[0,2].plot(0.10133, 0.206476, marker = 'X', c = 'grey', markersize = 20); ax[1,0].plot(0.10133, 0.206476, marker = 'X', c = 'grey', markersize = 20)
#ax[1,1].plot(0.10133, 0.206476, marker = 'X', c = 'grey', markersize = 20); ax[1,2].plot(0.10133, 0.206476, marker = 'X', c = 'grey', markersize = 20)
#plot rest wavelength with special marker:
ax[0,0].plot(e1['q'][16], e1['u'][16], marker = 'X', c = 'lime', markersize = 20, label = "rest\nwavelength")
ax[0,1].plot(e2['q'][16], e2['u'][16], marker = 'X', c = 'lime', markersize = 20, label = "rest λ")
ax[0,2].plot(e3['q'][17], e3['u'][17], marker = 'X', c = 'lime', markersize = 20, label = "rest λ")
ax[1,0].plot(e4['q'][16], e4['u'][16], marker = 'X', c = 'lime', markersize = 20)
ax[1,1].plot(e5['q'][17], e5['u'][17], marker = 'X', c = 'lime', markersize = 20)
ax[1,2].plot(e6['q'][17], e6['u'][17], marker = 'X', c = 'lime', markersize = 20)
#epoch 3&4 best fit minus outlyers
#e3bfit = get_lines(epoch3_20, 5530, 6210); e3bfit.remove_rows(slice(5, 9))
#e4bfit = get_lines(epoch4_20, 5530, 6210); e4bfit.remove_rows(slice(5, 7)) 
#bfit3 = best_fit(e3bfit['q'], e3bfit['u'], e3bfit['qerr'], e3bfit['uerr'])
#bfit4 = best_fit(e4bfit['q'], e4bfit['u'], e4bfit['qerr'], e4bfit['uerr'])
#QU data and label
epoch1 = QU(e1, ax[0,0]); epoch2 = QU(e2, ax[0,1]); epoch3 = QU(e3, ax[0,2])
epoch4 = QU(e4, ax[1,0]); epoch5 = QU(e5, ax[1,1]); epoch6 = QU(e6, ax[1,2])
cbar_ax = fig.add_axes([0.82, 0.15, 0.01, 0.7])
fig.colorbar(epoch1, cax=cbar_ax).set_label('(Å)', rotation= 0, fontsize = 20)
ax[0,0].text(1.6, 2.4, 'Days 0-7'); ax[0,1].text(1.6, 2.4, 'Day 26');ax[0,2].text(1.6, 2.4, 'Days 35-40')
ax[1,0].text(1.6, 2.4, 'Days 57-67'); ax[1,1].text(1.6, 2.4, 'Days 85-90');ax[1,2].text(1.6, 2.4, 'Day 295')
#plot circle where ISP could be to move orgin enough so that the feature doesn't cross the axis
ax[0,0].add_patch(plt.Circle((0, 0), radius = .32, edgecolor='m', facecolor = 'none', label = "ISP upper limit"))
ax[0,1].add_patch(plt.Circle((0, 0), radius = .32, edgecolor='m', facecolor = 'none'))
ax[0,2].add_patch(plt.Circle((0, 0), radius = .32, edgecolor='m', facecolor = 'none'))
ax[1,0].add_patch(plt.Circle((0, 0), radius = .32, edgecolor='m', facecolor = 'none'))
ax[1,1].add_patch(plt.Circle((0, 0), radius = .32, edgecolor='m', facecolor = 'none'))
ax[1,2].add_patch(plt.Circle((0, 0), radius = .32, edgecolor='m', facecolor = 'none'))
ax[0,0].legend(loc = 'lower left')
"""

#-------MYSTERY FEATURE (Epoch 6: 6500-7200)----------------------------------------
#plot_p_single([5876, 6600, 6678, 6731, 7065], ["He I", "Hα + N II", "He I", "S II", "He I"], ["b", "c", "b","r", "b"], data6, epoch6_20, 20, 10**-3.8, "Epoch6(Day295)")
#He_e6 = plot_line_velspace(epoch6_20, 5876, 15000, 10**-4.7, "Epoch6(Day295)")
#Mys_e6 = plot_line_velspace(epoch6_20, 7065, 15000, 10**-4.7, "Epoch6(Day295)", plot = 'on')

###Comparison of He 5876 & He 7065 (region blueward of rest)
He_reg = [[5572, 5876, 6180], [6700,7065,7430]]
He_reg_names = ["He 5876", "He7065"]
plot_QU_spec_comp2(unbindata = data6, bindata = epoch6_20, flx_adjust = .00015, region_lists = He_reg, region_names = He_reg_names, comp_color = "blue")

"""
plot_p_single([5876, 5762, 5701, 7065, 6853, 6927, 6051, 5989], ["He I", "-5816 & -8979km/s", "", "He I", "", "", "+5816 & +8979"], ["k", 'b', 'b', "k", "b", "b" , "r", "r"], data6, epoch6_20, 20, 10**-3.8, "Epoch6(Day295)")
Mys_e6blue = get_lines(epoch6_20, 6690, 7065)
Mys_e6blue_vel = get_vel_column(Mys_e6blue, 7065)
He5876_e6blue = get_lines(epoch6_20, 5570, 5900)
He5876_e6blue_vel = get_vel_column(He5876_e6blue, 5876)
"""
"""
He5876_e6red = get_lines(epoch6_20, 5870, 6180)
He5876_e6red_vel = get_vel_column(He5876_e6red, 5876)
#plt.step(He5876_e6blue['wave'], He5876_e6blue['p'], c = 'm', label = '-15600km/s from He 5876')
#plt.step(Mys_e6blue['wave'], Mys_e6blue['p'], c = 'g', label = '-15600km/s from He 7065')
#plt.legend()
#print(get_obswave(5876, -5816))
#print(get_velocity(5876, 5570));print(get_velocity(5876, 5700));print(get_velocity(5876, 5762)) #He 5876 line shifts
#print(get_velocity(7065,6697));print(get_obswave(7065, -8980)); print(get_obswave(7065, -5816)) #Wave values for velocities from He5876 shift applied to He 7065
#QU velocity space  
fig, ax = plt.subplots(1, figsize = (10, 10), sharex = True, sharey = True); fig.subplots_adjust(right = .85,  wspace=-.33, hspace=0.05)
fig.text(0.5, 0.08, '% Q Polarization', ha='center'); fig.text(0.15, 0.5, '% U Polarization', va='center', rotation='vertical') #set shared axis labels 
#epoch6 = QU(Mys_e6blue_vel, color_data = Mys_e6blue_vel['vel'],cmap_choice = 'autumn', bfit_c = 'm')
#epoch6_he = QU(He5876_e6blue_vel, color_data = He5876_e6blue_vel['vel'], cmap_choice = 'winter', bfit_c = 'b')
epoch6_he = QU(He5876_e6red_vel, color_data = He5876_e6red_vel['vel'], cmap_choice = 'spring', bfit_c = 'b')
#cbar_ax = fig.add_axes([0.9, 0.15, 0.01, 0.7])
#fig.colorbar(epoch6, cax=cbar_ax).set_label('Blueward of He 7065 (km/s)', rotation= 90, fontsize = 12)
cbar_ax = fig.add_axes([0.8, 0.15, 0.01, 0.7])
fig.colorbar(epoch6_he, cax=cbar_ax).set_label('redward of He 5876 (km/s)', rotation= 90, fontsize = 12)
ax.add_patch(plt.Circle((0, 0), radius = .32, edgecolor='red', facecolor = 'none', label = "ISP upper limit"))
"""
"""
fig, ax = plt.subplots(1, figsize = (10, 10), sharex = True, sharey = True); fig.subplots_adjust(right = .85,  wspace=-.33, hspace=0.05)
fig.text(0.5, 0.08, '% Q Polarization', ha='center'); fig.text(0.15, 0.5, '% U Polarization', va='center', rotation='vertical') #set shared axis labels 
epoch6 = QU(Mys_e6, cmap_choice = 'autumn', bfit_c = 'm')
epoch6_he = QU(He_e6, cmap_choice = 'winter', bfit_c = 'b')
cbar_ax = fig.add_axes([0.87, 0.15, 0.01, 0.7])
fig.colorbar(epoch6, cax=cbar_ax).set_label('(Å)', rotation= 0, fontsize = 20)
cbar_ax = fig.add_axes([0.8, 0.15, 0.01, 0.7])
fig.colorbar(epoch6_he, cax=cbar_ax).set_label('(Å)', rotation= 0, fontsize = 20)
ax.add_patch(plt.Circle((0, 0), radius = .32, edgecolor='red', facecolor = 'none', label = "ISP upper limit"))
"""
#-----------HE I 4471 & MG 4571 LINES------------------
"""
#all epochs 4471 & 4571 line region highlighted------------------------------------------------------------------------------------------
all_epochs20 = [epoch1_20, epoch2_20, epoch3_20, epoch4_20, epoch5_20, epoch6_20[0:-12]]
flx_adjust = [.00005, .00002, .0002, .0003, .0004, .00015]
epoch_names = ['days 0-7', 'day 26', 'days 35-40', 'days 57-67', 'days 85-90', 'day 295']
print(get_obswave(4571, 15000))
plot_all_epochs(all_epochs20, flx_adjust, epoch_names, xmin = 3700, regions_list = [[4342, 4799]], title = "SN 2012au He I & Mg +/-15,000km/s", line_list = [4471, 4571], color_list = ['b', 'm'], label_list = ['He I', 'Mg 4571'])
"""
"""
#Velocity Space plots - line Mg 4571 -------------------------------------------------------------------
e1 = plot_line_velspace(epoch1_20, 4571, 15000, 10**-4.5, "Epoch1(Day0-7)", save = 'on')
e2 = plot_line_velspace(epoch2_20, 4571, 15000, 10**-4.8, "Epoch2(Day26)", save = 'on')
e3 = plot_line_velspace(epoch3_20, 4571, 15000, 10**-4, "Epoch3(Day35-40)", save = 'on')
e4 = plot_line_velspace(epoch4_20, 4571, 15000, 10**-3.7, "Epoch4(Day57-67)", save = 'on')
e5 = plot_line_velspace(epoch5_20, 4571, 15000, 10**-3.5, "Epoch5(Day85-90)", save = 'on')
e6 = plot_line_velspace(epoch6_20, 4571, 15000, 10**-3.5, "Epoch6(Day295)", save = 'on')
"""
"""
#Single epoch highlight lines He I 4471 & Mg 4571-------------------------------------------------------------------------------------
Mg = [4471, 4571]; Mg_lab = ["He I (4471)", "Mg (4571)"]; Mg_c = ["b", "m"]
#epoch1
plot_p_single(Mg, Mg_lab, Mg_c, unbin1, epoch1_20, 20, 10**-4.5, "Epoch_1")
plt.plot(e1['wave'], e1['p'], c= 'y', label = "+/- 20000km/s");plt.legend(loc="upper right", fontsize = 12)
plt.savefig("Mg_4571/Epoch1_binned20")
#epoch2
plot_p_single(Mg, Mg_lab, Mg_c, data2, epoch2_20, 20, 10**-4.7, "Epoch_2")
plt.plot(e2['wave'], e2['p'], c= 'y', label = "+/- 20000km/s");plt.legend(loc="upper right", fontsize = 12)
plt.savefig("Mg_4571/Epoch2_binned20")
#epoch3
plot_p_single(Mg, Mg_lab, Mg_c, data3, epoch3_20, 20, 10**-3.6, "Epoch_3")
plt.plot(e3['wave'], e3['p'], c= 'y', label = "+/- 20000km/s");plt.legend(loc="upper right", fontsize = 12)
plt.savefig("Mg_4571/Epoch3_binned20")
#epoch4
plot_p_single(Mg, Mg_lab, Mg_c, data4, epoch4_20, 20, 10**-3.5, "Epoch_4")
plt.plot(e4['wave'], e4['p'],c= 'y', label = "+/- 20000km/s");plt.legend(loc="upper right", fontsize = 12)
plt.savefig("Mg_4571/Epoch4_binned20")
#epoch5
plot_p_single(Mg, Mg_lab, Mg_c, data5, epoch5_20, 20, 10**-3.2, "Epoch_5")
plt.plot(e5['wave'], e5['p'],c= 'y', label = "+/- 20000km/s");plt.legend(loc="upper right", fontsize = 12)
plt.savefig("Mg_4571/Epoch5_binned20")
#epoch6
plot_p_single(Mg, Mg_lab, Mg_c, data6, epoch6_20, 20, 10**-3.8, "Epoch_6")
plt.plot(e6['wave'], e6['p'], c= 'y', label = "+/- 20000km/s"); plt.legend(loc="upper right", fontsize = 12)
plt.savefig("Mg_4571/Epoch6_binned20")
"""
"""
#QU all epoch subplots combined - line He 4471/ Mg 4571--------------------------------------------------------------------------------------
import matplotlib.pyplot as plt
fig, ax = plt.subplots(2, 3, figsize = (30, 16)); fig.subplots_adjust(right=0.8)
#local continuum region data
all_epochs = [epoch1_20, epoch2_20, epoch3_20, epoch4_20, epoch5_20, epoch6_20[0:-12]]
cont1 = get_aves(all_epochs, [5100, 5500]) 
qcont = cont1['q']; ucont = cont1['u']; qcerr = cont1['qerr']; ucerr = cont1['uerr']
ax[0,0].plot(qcont[0], ucont[0], marker = '*', c = 'k', markersize = 20); ax[0,1].plot(qcont[1], ucont[1], marker = '*', c = 'k', markersize = 20)
ax[0,2].plot(qcont[2], ucont[2], marker = '*', c = 'k', markersize = 20); ax[1,0].plot(qcont[3], ucont[3], marker = '*', c = 'k', markersize = 20)
ax[1,1].plot(qcont[4], ucont[4], marker = '*', c = 'k', markersize = 20); ax[1,2].plot(qcont[5], ucont[5], marker = '*', c = 'k', markersize = 20)
#epoch 3&4 best fit minus outlyers
#e3bfit = get_lines(epoch3_20, 5530, 6210); e3bfit.remove_rows(slice(5, 9))
#e4bfit = get_lines(epoch4_20, 5530, 6210); e4bfit.remove_rows(slice(5, 7)) 
#bfit3 = best_fit(e3bfit['q'], e3bfit['u'], e3bfit['qerr'], e3bfit['uerr'])
#bfit4 = best_fit(e4bfit['q'], e4bfit['u'], e4bfit['qerr'], e4bfit['uerr'])
#QU data and label
epoch1 = QU(e1, ax[0,0]); epoch2 = QU(e2, ax[0,1]); epoch3 = QU(e3, ax[0,2])
epoch4 = QU(e4, ax[1,0]); epoch5 = QU(e5, ax[1,1]); epoch6 = QU(e6, ax[1,2])
cbar_ax = fig.add_axes([0.82, 0.15, 0.01, 0.7])
fig.colorbar(epoch1, cax=cbar_ax).set_label('(Å)', rotation= 0, fontsize = 20)
ax[0,0].text(0.016, 0.034, 'Days 0-7'); ax[0,1].text(0.016, 0.034, 'Day 26');ax[0,2].text(0.016, 0.034, 'Days 35-40')
ax[1,0].text(0.016, 0.034, 'Days 57-67'); ax[1,1].text(0.016, 0.034, 'Days 85-90');ax[1,2].text(0.016, 0.034, 'Day 295')
"""
#-----------Fe II 4924 LINE------------------
"""
#all epochs 4924 line region highlighted------------------------------------------------------------------------------------------
all_epochs20 = [epoch1_20, epoch2_20, epoch3_20, epoch4_20, epoch5_20, epoch6_20[0:-12]]
flx_adjust = [.00005, .00002, .0002, .0003, .0004, .00015]
epoch_names = ['days 0-7', 'day 26', 'days 35-40', 'days 57-67', 'days 85-90', 'day 295']
#print(get_obswave(4924, 20000))
plot_all_epochs(all_epochs20, flx_adjust, epoch_names, xmin = 3700, regions_list = [[4595, 5252]], title = "SN 2012au Fe II +/- 20,000km/s", line_list = [4924], color_list = ['g'], label_list = ['Fe II (4924)'])
"""
"""
#plot of epochs with Fe II velocity shift labeled ------------------------------------------------------------
epochs1_4_20 = [epoch1_20, epoch2_20, epoch3_20, epoch4_20, epoch5_20, epoch6_20[0:-12]]
flx_adjust = [.00005, .00002, .0002, .0003, .0004, .00015]
epoch_names = ['days 0-7', 'day 26', 'days 35-40', 'days 57-67', 'days 85-90', 'day 295']
im = plot_all_epochs(epochs1_4_20, flx_adjust, epoch_names, xmin = 3700, line_list=[4924], color_list= ['g'], label_list= ["Fe II 4924"]) 
im.axvline(x = 5582, ymin = .70, ymax = .98, c = 'c'); im.text(5250, 16.5, '-15,000 km/s') #velocities based on flux minimums
im.axvline(x = 5650, ymin = .35, ymax = .65, c = 'c'); im.text(5250, 11, '-11,530 km/s')
im.axvline(x = 5665, ymin = .05, ymax = .30, c = 'c'); im.text(5250, 5.5, '-10,765 km/s')
"""
"""
#Velocity Space plots - line Fe II 4924 -------------------------------------------------------------------
e1 = plot_line_velspace(epoch1_20, 4924, 20000, 10**-4.5, "Epoch1(Day0-7)");
e2 = plot_line_velspace(epoch2_20, 4924, 20000, 10**-4.8, "Epoch2(Day26)"); 
e3 = plot_line_velspace(epoch3_20, 4924, 20000, 10**-4, "Epoch3(Day35-40)")
e4 = plot_line_velspace(epoch4_20, 4924, 20000, 10**-3.7, "Epoch4(Day57-67)")
e5 = plot_line_velspace(epoch5_20, 4924, 20000, 10**-3.5, "Epoch5(Day85-90)")
e6 = plot_line_velspace(epoch6_20, 4924, 20000, 10**-3.5, "Epoch6(Day295)")
"""
"""
#stacked flux and polarization
epochs = [e1, e2, e3, e4, e5, e6]
flx_adjust= [.00005, .00002, .0004, .0007, .0012, .0001]
epoch_names = ['days 0-7', 'day 26', 'days 35-40', 'days 57-67', 'days 85-90', 'day 295']

plot_all_epochs(epochs, flx_adjust, epoch_names, title= 'Fe II Rest λ (4924)', line_list= [4924], label_list= [''], color_list= ['k'])
"""
"""
#Single epoch highlight line Fe II 4924-------------------------------------------------------------------------------------
Fe = [4924]; Fe_lab = ["Fe II (4924)"]; Fe_c = ["g"]
#epoch1
plot_p_single(Fe, Fe_lab, Fe_c, unbin1, epoch1_20, 20, 10**-4.5, "Epoch_1")
plt.plot(e1['wave'], e1['p'], c= 'b', label = "+/- 20000km/s");plt.legend(loc="upper right", fontsize = 12)
plt.savefig("FeII_4924/Epoch1_binned20")
#epoch2
plot_p_single(Fe, Fe_lab, Fe_c, data2, epoch2_20, 20, 10**-4.7, "Epoch_2")
plt.plot(e2['wave'], e2['p'], c= 'b', label = "+/- 20000km/s");plt.legend(loc="upper right", fontsize = 12)
plt.savefig("FeII_4924/Epoch2_binned20")
#epoch3
plot_p_single(Fe, Fe_lab, Fe_c, data3, epoch3_20, 20, 10**-3.6, "Epoch_3")
plt.plot(e3['wave'], e3['p'], c= 'b', label = "+/- 20000km/s");plt.legend(loc="upper right", fontsize = 12)
plt.savefig("FeII_4924/Epoch3_binned20")
#epoch4
plot_p_single(Fe, Fe_lab, Fe_c, data4, epoch4_20, 20, 10**-3.5, "Epoch_4")
plt.plot(e4['wave'], e4['p'],c= 'b', label = "+/- 20000km/s");plt.legend(loc="upper right", fontsize = 12)
plt.savefig("FeII_4924/Epoch4_binned20")
#epoch5
plot_p_single(Fe, Fe_lab, Fe_c, data5, epoch5_20, 20, 10**-3.2, "Epoch_5")
plt.plot(e5['wave'], e5['p'],c= 'b', label = "+/- 20000km/s");plt.legend(loc="upper right", fontsize = 12)
plt.savefig("FeII_4924/Epoch5_binned20")
#epoch6
plot_p_single(Fe, Fe_lab, Fe_c, data6, epoch6_20, 20, 10**-3.8, "Epoch_6")
plt.plot(e6['wave'], e6['p'], c= 'b', label = "+/- 20000km/s"); plt.legend(loc="upper right", fontsize = 12)
plt.savefig("FeII_4924/Epoch6_binned20")
"""
"""
#QU all epoch subplots combined - line Fe II 4924--------------------------------------------------------------------------------------
import matplotlib.pyplot as plt
fig, ax = plt.subplots(2, 3, figsize = (20, 13)); fig.subplots_adjust(right=0.8)
#local continuum region data
all_epochs = [epoch1_20, epoch2_20, epoch3_20, epoch4_20, epoch5_20, epoch6_20[0:-12]]
cont1 = get_aves(all_epochs, [5100, 5500]) 
qcont = cont1['q']; ucont = cont1['u']; qcerr = cont1['qerr']; ucerr = cont1['uerr']
ax[0,0].plot(qcont[0], ucont[0], marker = '*', c = 'k', markersize = 20); ax[0,1].plot(qcont[1], ucont[1], marker = '*', c = 'k', markersize = 20)
ax[0,2].plot(qcont[2], ucont[2], marker = '*', c = 'k', markersize = 20); ax[1,0].plot(qcont[3], ucont[3], marker = '*', c = 'k', markersize = 20)
ax[1,1].plot(qcont[4], ucont[4], marker = '*', c = 'k', markersize = 20); ax[1,2].plot(qcont[5], ucont[5], marker = '*', c = 'k', markersize = 20)
#Pandey et al. 21 ISP estimate as grey cross
ax[0,0].plot(0.10133, 0.206476, marker = 'X', c = 'grey', markersize = 20); ax[0,1].plot(0.10133, 0.206476, marker = 'X', c = 'grey', markersize = 20)
ax[0,2].plot(0.10133, 0.206476, marker = 'X', c = 'grey', markersize = 20); ax[1,0].plot(0.10133, 0.206476, marker = 'X', c = 'grey', markersize = 20)
ax[1,1].plot(0.10133, 0.206476, marker = 'X', c = 'grey', markersize = 20); ax[1,2].plot(0.10133, 0.206476, marker = 'X', c = 'grey', markersize = 20)
#plot rest wavelength with special marker:
ax[0,0].plot(e1['q'][16], e1['u'][16], marker = 'X', c = 'lime', markersize = 20, label = "rest λ")
ax[0,1].plot(e2['q'][16], e2['u'][16], marker = 'X', c = 'lime', markersize = 20)
ax[0,2].plot(e3['q'][16], e3['u'][16], marker = 'X', c = 'lime', markersize = 20)
ax[1,0].plot(e4['q'][16], e4['u'][16], marker = 'X', c = 'lime', markersize = 20)
ax[1,1].plot(e5['q'][16], e5['u'][16], marker = 'X', c = 'lime', markersize = 20)
ax[1,2].plot(e6['q'][16], e6['u'][16], marker = 'X', c = 'lime', markersize = 20)
#plot our ISP estimate from epoch 5 continuum regions weighted average:
qisp = 0.1716573487345389; uisp = 0.18495597891967216
ax[0,0].plot(qisp, uisp, marker = 'X', c = 'k', markersize = 20)
ax[0,1].plot(qisp, uisp, marker = 'X', c = 'k', markersize = 20)
ax[0,2].plot(qisp, uisp, marker = 'X', c = 'k', markersize = 20)
ax[1,0].plot(qisp, uisp, marker = 'X', c = 'k', markersize = 20)
ax[1,1].plot(qisp, uisp, marker = 'X', c = 'k', markersize = 20)
ax[1,2].plot(qisp, uisp, marker = 'X', c = 'k', markersize = 20)

#epoch 3&4 best fit minus outlyers
#e3bfit = get_lines(epoch3_20, 5530, 6210); e3bfit.remove_rows(slice(5, 9))
#e4bfit = get_lines(epoch4_20, 5530, 6210); e4bfit.remove_rows(slice(5, 7)) 
#bfit3 = best_fit(e3bfit['q'], e3bfit['u'], e3bfit['qerr'], e3bfit['uerr'])
#bfit4 = best_fit(e4bfit['q'], e4bfit['u'], e4bfit['qerr'], e4bfit['uerr'])
#QU data and label
epoch1 = QU(e1, ax[0,0]); epoch2 = QU(e2, ax[0,1]); epoch3 = QU(e3, ax[0,2])
epoch4 = QU(e4, ax[1,0]); epoch5 = QU(e5, ax[1,1]); epoch6 = QU(e6, ax[1,2])
cbar_ax = fig.add_axes([0.82, 0.15, 0.01, 0.7])
fig.colorbar(epoch1, cax=cbar_ax).set_label('(Å)', rotation= 0, fontsize = 16)
ax[0,0].text(-2.4, 2.4, 'Days 0-7'); ax[0,1].text(-2.4, 2.4, 'Day 26');ax[0,2].text(-2.4, 2.4, 'Days 35-40')
ax[1,0].text(-2.4, 2.4, 'Days 57-67'); ax[1,1].text(-2.4, 2.4, 'Days 85-90');ax[1,2].text(-2.4, 2.4, 'Day 295')
"""

#-------Fe II 5169 LINE----------------------------------------
"""
#all epochs 5169 line region highlighted------------------------------------------------------------------------------------------
all_epochs20 = [epoch1_20, epoch2_20, epoch3_20, epoch4_20, epoch5_20, epoch6_20[0:-12]]
flx_adjust = [.00005, .00002, .0002, .0003, .0004, .00015]
epoch_names = ['days 0-7', 'day 26', 'days 35-40', 'days 57-67', 'days 85-90', 'day 295']
#print(get_obswave(5169, 20000))
plot_all_epochs(all_epochs20, flx_adjust, epoch_names, xmin = 3700, regions_list = [[4824, 5514]], title = "SN 2012au Fe II +/- 20,000km/s", line_list = [5169], color_list = ['g'], label_list = ['Fe II (5169)'])
"""
"""
#Velocity Space plots - line Fe II 5169 -------------------------------------------------------------------
e1 = plot_line_velspace(epoch1_20, 5169, 20000, 10**-4.5, "Epoch1(Day0-7)");
e2 = plot_line_velspace(epoch2_20, 5169, 20000, 10**-4.8, "Epoch2(Day26)"); 
e3 = plot_line_velspace(epoch3_20, 5169, 20000, 10**-4, "Epoch3(Day35-40)")
e4 = plot_line_velspace(epoch4_20, 5169, 20000, 10**-3.7, "Epoch4(Day57-67)")
e5 = plot_line_velspace(epoch5_20, 5169, 20000, 10**-3.5, "Epoch5(Day85-90)")
e6 = plot_line_velspace(epoch6_20, 5169, 20000, 10**-3.5, "Epoch6(Day295)")
"""
"""
#Single epoch highlight line Fe II 5169-------------------------------------------------------------------------------------
Fe = [5169]; Fe_lab = ["Fe II (5169)"]; Fe_c = ["g"]
#epoch1
plot_p_single(Fe, Fe_lab, Fe_c, unbin1, epoch1_20, 20, 10**-4.5, "Epoch_1")
plt.plot(e1['wave'], e1['p'], c= 'b', label = "+/- 20000km/s");plt.legend(loc="upper right", fontsize = 12)
plt.savefig("FeII_5169/Epoch1_binned20")
#epoch2
plot_p_single(Fe, Fe_lab, Fe_c, data2, epoch2_20, 20, 10**-4.7, "Epoch_2")
plt.plot(e2['wave'], e2['p'], c= 'b', label = "+/- 20000km/s");plt.legend(loc="upper right", fontsize = 12)
plt.savefig("FeII_5169/Epoch2_binned20")
#epoch3
plot_p_single(Fe, Fe_lab, Fe_c, data3, epoch3_20, 20, 10**-3.6, "Epoch_3")
plt.plot(e3['wave'], e3['p'], c= 'b', label = "+/- 20000km/s");plt.legend(loc="upper right", fontsize = 12)
plt.savefig("FeII_5169/Epoch3_binned20")
#epoch4
plot_p_single(Fe, Fe_lab, Fe_c, data4, epoch4_20, 20, 10**-3.5, "Epoch_4")
plt.plot(e4['wave'], e4['p'],c= 'b', label = "+/- 20000km/s");plt.legend(loc="upper right", fontsize = 12)
plt.savefig("FeII_5169/Epoch4_binned20")
#epoch5
plot_p_single(Fe, Fe_lab, Fe_c, data5, epoch5_20, 20, 10**-3.2, "Epoch_5")
plt.plot(e5['wave'], e5['p'],c= 'b', label = "+/- 20000km/s");plt.legend(loc="upper right", fontsize = 12)
plt.savefig("FeII_5169/Epoch5_binned20")
#epoch6
plot_p_single(Fe, Fe_lab, Fe_c, data6, epoch6_20, 20, 10**-3.8, "Epoch_6")
plt.plot(e6['wave'], e6['p'], c= 'b', label = "+/- 20000km/s"); plt.legend(loc="upper right", fontsize = 12)
plt.savefig("FeII_5169/Epoch6_binned20")
"""
"""
#QU all epoch subplots combined - line Fe II 5169--------------------------------------------------------------------------------------
import matplotlib.pyplot as plt
fig, ax = plt.subplots(2, 3, figsize = (30, 16)); fig.subplots_adjust(right=0.8)
#local continuum region data
all_epochs = [epoch1_20, epoch2_20, epoch3_20, epoch4_20, epoch5_20, epoch6_20[0:-12]]
cont1 = get_aves(all_epochs, [5100, 5500]) 
qcont = cont1['q']; ucont = cont1['u']; qcerr = cont1['qerr']; ucerr = cont1['uerr']
ax[0,0].plot(qcont[0], ucont[0], marker = '*', c = 'k', markersize = 20); ax[0,1].plot(qcont[1], ucont[1], marker = '*', c = 'k', markersize = 20)
ax[0,2].plot(qcont[2], ucont[2], marker = '*', c = 'k', markersize = 20); ax[1,0].plot(qcont[3], ucont[3], marker = '*', c = 'k', markersize = 20)
ax[1,1].plot(qcont[4], ucont[4], marker = '*', c = 'k', markersize = 20); ax[1,2].plot(qcont[5], ucont[5], marker = '*', c = 'k', markersize = 20)
#epoch 3&4 best fit minus outlyers
#e3bfit = get_lines(epoch3_20, 5530, 6210); e3bfit.remove_rows(slice(5, 9))
#e4bfit = get_lines(epoch4_20, 5530, 6210); e4bfit.remove_rows(slice(5, 7)) 
#bfit3 = best_fit(e3bfit['q'], e3bfit['u'], e3bfit['qerr'], e3bfit['uerr'])
#bfit4 = best_fit(e4bfit['q'], e4bfit['u'], e4bfit['qerr'], e4bfit['uerr'])
#QU data and label
epoch1 = QU(e1, ax[0,0]); epoch2 = QU(e2, ax[0,1]); epoch3 = QU(e3, ax[0,2])
epoch4 = QU(e4, ax[1,0]); epoch5 = QU(e5, ax[1,1]); epoch6 = QU(e6, ax[1,2])
cbar_ax = fig.add_axes([0.82, 0.15, 0.01, 0.7])
fig.colorbar(epoch1, cax=cbar_ax).set_label('(Å)', rotation= 0, fontsize = 20)
ax[0,0].text(0.016, 0.034, 'Days 0-7'); ax[0,1].text(0.016, 0.034, 'Day 26');ax[0,2].text(0.016, 0.034, 'Days 35-40')
ax[1,0].text(0.016, 0.034, 'Days 57-67'); ax[1,1].text(0.016, 0.034, 'Days 85-90');ax[1,2].text(0.016, 0.034, 'Day 295')
"""

#-----------Fe II 4924 & He I 5876 LINE------------------
"""
#all epochs 4924 line region highlighted------------------------------------------------------------------------------------------
all_epochs20 = [epoch1_20, epoch2_20, epoch3_20, epoch4_20, epoch5_20, epoch6_20[0:-12]]
flx_adjust = [.00005, .00002, .0002, .0003, .0004, .00015]
epoch_names = ['days 0-7', 'day 26', 'days 35-40', 'days 57-67', 'days 85-90', 'day 295']
im = plot_all_epochs(all_epochs20, flx_adjust, epoch_names, xmin = 3700, regions_list = [[4678, 5170], [5582, 6170]], line_list = [4924, 5876], color_list = ['g', 'b'], label_list = ['Fe II (4924)', 'He I (5876)'])
im.axvline(x = 5582, ymin = .85, ymax = .98, c = 'k'); im.text(5250, 30, '-15,000 km/s') #velocities based on flux minimums
im.axvline(x = 5650, ymin = .69, ymax = .82, c = 'k'); im.text(5250, 25, '-11,530 km/s')
im.axvline(x = 5665, ymin = .52, ymax = .65, c = 'k'); im.text(5250, 20, '-10,765 km/s')
im.axvline(x = 4705, ymin = .85, ymax = .98, c = 'k', linestyle = "--"); im.text(4350, 30, '-13,333 km/s')
im.axvline(x = 4750, ymin = .69, ymax = .82, c = 'k', linestyle = "--"); im.text(4350, 25, '-10,594 km/s')
im.axvline(x = 4760, ymin = .52, ymax = .65, c = 'k', linestyle = "--"); im.text(4370, 20, '-9,984 km/s')
"""
"""
#Velocity space plots with shifted flux features labeled ------------------------------------------------------------------------
#epoch1
line = [4678, 4924, 5582, 5876]; lab = ["","Fe II (4924)", "", "He I (5876)"]; c = ["springgreen","g", "royalblue", "b"]
plot_p_single(line, lab, c, unbin1, epoch1_20, 20, 10**-4.5, "Epoch_1 (Shifted -15,000km/s)")
#epoch2
line = [4735, 4924, 5650, 5876]; lab = ["","Fe II (4924)", "", "He 2 (5876)"]; c = ["springgreen","g", "royalblue", "b"] 
plot_p_single(line, lab, c, data2, epoch2_20, 20, 10**-4.7, "Epoch_2 (Shifted -11,530km/s)")
#epoch3
line = [4747, 4924, 5665, 5876]; lab = ["","Fe II (4924)", "", "He 2 (5876)"]; c = ["springgreen","g", "royalblue", "b"] 
plot_p_single(line, lab, c, data3, epoch3_20, 20, 10**-3.6, "Epoch_3 (Shifted -10,765km/s)")
"""
"""
#unfinished***** code to stack PA, P and FLX for two lines being compared
rawdat = [unbin1, data2, data3, data4, data5]
count = 1
fig, (ax0, ax1, ax2) = plt.subplots(3, figsize = (10, 14), sharex = True)
for data in rawdat:
    epoch = Bin_data(data, 4, 20) 
    Fe = get_lines(epoch, 4595, 5252); Fe_raw = get_lines(data, 4595, 5252)
    He = get_lines(epoch, 5484, 6268); He_raw = get_lines(data, 5484, 6268)
    Fe_vels = []; Fe_vels_raw = []; He_vels = []; He_vels_raw = []
    Fe_PA = PA(Fe['q'], Fe['u']); He_PA = PA(He['q'], He['u'])
    for i in Fe['wave']:
        Fe_vels.append(get_velocity(4924, i))
    for i in Fe_raw['wave']:
        Fe_vels_raw.append(get_velocity(4924, i))
    for i in He['wave']:
        He_vels.append(get_velocity(5976, i))
    for i in He_raw['wave']:
        He_vels_raw.append(get_velocity(5976, i))
   
    #fig, (ax0, ax1, ax2) = plt.subplots(3, figsize = (5, 7), sharex = True)
    ax0.set_ylabel('Rel. Flux'); ax0.plot(Fe_vels_raw, Fe_raw['flx'], label = 'Fe 4924 epoch' + str(count), linestyle='dashed'); ax0.plot(He_vels_raw, He_raw['flx'], label = 'He 5876')
    ax1.set_ylabel('Total Polarization (%)'); ax1.step(Fe_vels, Fe['p'], linestyle='dashed'); ax1.step(He_vels, He['p'])
    ax2.set_ylabel('PA (°)'); ax2.set_xlabel('Velocity (km/s)'); ax2.step(Fe_vels, Fe_PA, linestyle='dashed'); ax2.step(He_vels, He_PA)
    ax0.axvline(0, c = 'grey'); ax1.axvline(0, c = 'grey'); ax2.axvline(0, c = 'grey');
    plt.subplots_adjust(hspace=0); plt.suptitle('Epoch '+str(count)); fig.legend(loc = 'upper right')
    
    count = count+1   
"""
#-----------All HE LINES------------------
#Single epoch highlight He lines-------------------------------------------------------------------------------------
"""
#epoch1 -flux shift
He = [4471, 5618, 5876, 6385, 6678, 6750, 7065 ]; He_lab = ["4471", "-13163km/s", "5876", "-13153km/s", "6678", "-13366km/s", "7065"]; He_c = ["k","b", "k", "b", "k","b", "k"]
plot_p_single(He, He_lab, He_c, unbin1, epoch1_20, 20, 10**-4.5, "Days 0-7 He I (flux shifts)")
plt.arrow(5876, 3, -200, 0, head_width = 0.08, head_length = 50)
plt.arrow(6678, 3.4, -220, 0, head_width = 0.08, head_length = 50)
plt.arrow(7065, 3.7, -250, 0, head_width = 0.08, head_length = 50)
maxs, mins = find_lines(unbin1['wave'], unbin1['flx'], 10, 50) #use find lines to get mins 
mins = [x[0] for x in mins]; c = ["b" for x in range(len(mins))]
plot_p_single(mins, mins, c,  unbin1, epoch1_20, 20, 10**-4.5, "Epoch_1") #plot find_lines mins to get values 
print(get_velocity(5876, 5618)) #calculate velocity shift from rest to mins to plot/label 
print(get_velocity(6678, 6385))
print(get_velocity(7065, 6750)) #got this shift from Pandey et al. 
#epoch1 -polarization shift
He_p = [4471, 5675, 5876, 6430, 6678, 6810, 7065 ]; He_p_lab = ["4471", "-10255km/s", "5876", "-11133km/s", "6678", "-10820km/s", "7065"]; He_p_c = ["k","b", "k", "b", "k","b", "k"]
plot_p_single(He_p, He_p_lab, He_p_c, unbin1, epoch1_20, 20, 10**-4.5, "Days 0-7 He I (polarization shifts)")
plt.arrow(5876, 3, -150, 0, head_width = 0.08, head_length = 50)
plt.arrow(6678, 3.4, -200, 0, head_width = 0.08, head_length = 50)
plt.arrow(7065, 3.7, -200, 0, head_width = 0.08, head_length = 50)
print(get_velocity(5876, 5675)) #calculate velocity shift from rest to mins to plot/label 
print(get_velocity(6678, 6430))
print(get_velocity(7065, 6810)) 
"""
"""#epoch2 - flux shift
He = [4471, 5675, 5876, 6460, 6678, 6865, 7065 ]; He_lab = ["4471", "-10255km/s", "5876", "-9787km/s", "6678", "-8487km/s", "7065"]; He_c = ["k","b", "k", "b", "k","b", "k"]
plot_p_single(He, He_lab, He_c, data2, epoch2_20, 20, 10**-4.8, "Day 26 He I (flux shift)")
maxs, mins = find_lines(data2['wave'], data2['flx'], 10, 50)
mins = [x[0] for x in mins]; c = ["b" for x in range(len(mins))]
#plot_p_single(mins, mins, c,  data2, epoch2_20, 20, 10**-4.8, "Epoch_2")
plt.arrow(5876, 3, -150, 0, head_width = 0.08, head_length = 50)
plt.arrow(6678, 3.4, -150, 0, head_width = 0.08, head_length = 50)
plt.arrow(7065, 3.7, -120, 0, head_width = 0.08, head_length = 50)
print(get_velocity(5876, 5675))
print(get_velocity(6678, 6460))#got this from Pandey et al. 
print(get_velocity(7065, 6865))#got this shift from Pandey et al.
"""
"""#epoch2 - Polarization shift
He = [4471, 5720, 5876, 6550, 6678, 6865, 7065 ]; He_lab = ["4471", "-7959km/s", "5876", "-5746km/s", "6678", "-8487km/s", "7065"]; He_c = ["k","b", "k", "b", "k","b", "k"]
plot_p_single(He, He_lab, He_c, data2, epoch2_20, 20, 10**-4.8, "Day 26 He I (Polarization shift)")
maxs, mins = find_lines(data2['wave'], data2['flx'], 10, 50)
mins = [x[0] for x in mins]; c = ["b" for x in range(len(mins))]
#plot_p_single(mins, mins, c,  data2, epoch2_20, 20, 10**-4.8, "Epoch_2")
plt.arrow(5876, 3, -100, 0, head_width = 0.08, head_length = 50)
plt.arrow(6678, 3.4, -80, 0, head_width = 0.08, head_length = 50)
plt.arrow(7065, 3.7, -120, 0, head_width = 0.08, head_length = 50)
print(get_velocity(5876, 5720))
print(get_velocity(6678, 6550)) 
print(get_velocity(7065, 6865))
"""
"""#epoch3 -flux shift 
He = [4471, 5675, 5876, 6470, 6678, 6880, 7065 ]; He_lab = ["4471", "-9337km/s", "5876", "-9338km/s", "6678","-7850km/s", "7065"]; He_c = ["k","b", "k", "b", "k","b", "k"]
plot_p_single(He, He_lab, He_c, data3, epoch3_20, 20, 10**-3.8, "Days 35-40 He I (flux shift)")
plt.arrow(5876, 3, -120, 0, head_width = 0.08, head_length = 50)
plt.arrow(6678, 3.4, -120, 0, head_width = 0.08, head_length = 50)
plt.arrow(7065, 3.7, -120, 0, head_width = 0.08, head_length = 50)
maxs, mins = find_lines(data3['wave'], data3['flx'], 10, 50)
mins = [x[0] for x in mins]; c = ["b" for x in range(len(mins))]
#plot_p_single(mins, mins, c,  data3, epoch3_20, 20, 10**-3.8, "Epoch_3")
print(get_velocity(5876, 5693))
print(get_velocity(6678, 6470)) #from Pandey et al. 
print(get_velocity(7065, 6880))
#epoch3 - polarization shift (did not see distinct features to label so no plot was created)
"""
"""#epoch4
He = [4471, 5707, 5876, 6510, 6678, 6880, 7065 ]; He_lab = ["4471", "-8622km/s", "5876", "-7541km/s", "6678","-7426km/s",  "7065"]; He_c = ["k","b", "k", "b", "k", "b", "k"]
plot_p_single(He, He_lab, He_c, data4, epoch4_20, 20, 10**-3.5, "Days 57-67 He I (flux shift)")
plt.arrow(5876, 3, -120, 0, head_width = 0.08, head_length = 50)
plt.arrow(6678, 3.4, -120, 0, head_width = 0.08, head_length = 50)
plt.arrow(7065, 3.7, -120, 0, head_width = 0.08, head_length = 50)
maxs, mins = find_lines(data4['wave'], data4['flx'], 10, 50)
mins = [x[0] for x in mins]; c = ["b" for x in range(len(mins))]
#plot_p_single(mins, mins, c,  data4, epoch4_20, 20, 10**-3.5, "Epoch_2")
print(get_velocity(5876, 5707))
print(get_velocity(6678, 6510))
print(get_velocity(7065, 6880))
#epoch4 - polarization shifts (did not see distinct features to label so no plot was created)
"""
"""#epoch5
He = [4471, 5711, 5876, 6510, 6678, 6885, 7065 ]; He_lab = ["4471", "-8418km/s", "5876", "-7542km/s", "6678", "-7638km/s", "7065"]; He_c = ["k","b", "k", "b", "k", "b", "k"]
plot_p_single(He, He_lab, He_c, data5, epoch5_20, 20, 10**-3.2, "Days 85-90 He I (flux shift)")
plt.arrow(5876, 3, -120, 0, head_width = 0.08, head_length = 50)
plt.arrow(6678, 3.4, -120, 0, head_width = 0.08, head_length = 50)
plt.arrow(7065, 3.7, -120, 0, head_width = 0.08, head_length = 50)
maxs, mins = find_lines(data5['wave'], data5['flx'], 10, 50)
mins = [x[0] for x in mins]; c = ["b" for x in range(len(mins))]
#plot_p_single(mins, mins, c,  data5, epoch5_20, 20, 10**-3.2, "Epoch_2")
print(get_velocity(5876, 5711))
print(get_velocity(6678, 6510))#from Pandey et al. 
print(get_velocity(7065, 6885))
#epoch4 - polarization shifts (did not see distinct features to label so no plot was created)
"""
"""#epoch6
He = [4471, 5810, 5930, 5876, 6678, 7030, 7065 ]; He_lab = ["4471", "", "-3367km/s \n 5876 \n +2755km/s", "", "6678", " ", "-1485km/s \n 7065"]; He_c = ["k","b","b", "k", "k", "b", "k"]
plot_p_single(He, He_lab, He_c, data6, epoch6_20, 20, 10**-3.8, "Day 295 He I (flux shift)")

print(get_velocity(5876, 5930)) #from Pandey et al. noted as +NaID 
print(get_velocity(6678, 6600))
print(get_velocity(5876, 5810))
"""
"""
#All epoch feature velocities
flux_features = [-13163, -10204, -9337, -8622, -3673, -3520]; pol_features = [-10051, -7143, 3367, 7347, -15969, 306]
pol_flux_features = []; days = [0, 26, 35, 57, 85, 295]
fig, (ax0, ax1, ax2) = plt.subplots(3, sharex = True, sharey = True, figsize = (10, 12)); 
ax0.scatter(days, flux_features); ax0.set_ylabel("FLux Feature Velocity (km/s)"); ax0.axhline(0, c = 'k', lw = 1)
ax1.scatter(days, pol_features); ax1.set_ylabel("Polarization Feature Velocity (km/s)"); ax1.axhline(0, c = 'k', lw = 1)
ax2.set_ylabel("Polarized Flux Feature Velocity (km/s)"); ax2.set_xlabel("Days Since Max Brightness")
"""