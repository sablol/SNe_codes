# -*- coding: utf-8 -*-
"""
Created on Thu Aug 26 10:28:59 2021

@author: sabri
"""
from get_data import * 
#from spec_bin import * 
from Binning_routines import * 
import matplotlib.pyplot as plt
from astropy.io import fits as pyfits
"""
#unbin = get_txtFITS("SN_2010jl", "sn2010jl_1_allcomb.flx.txt", "sn2010jl_1_allcomb.q.txt", "sn2010jl_1_allcomb.qsig.txt", "sn2010jl_1_allcomb.flx.txt", "sn2010jl_1_allcomb.u.txt", "sn2010jl_1_allcomb.usig.txt")
unbin = get_txtFITS("Epoch_1/all_comb", "sn2012au_allcomb.flx.txt" , "sn2012au_allcomb.q.txt", "sn2012au_allcomb.qsig.txt", "sn2012au_allcomb.qsum.txt", "sn2012au_allcomb.u.txt", "sn2012au_allcomb.usig.txt")
#unbin = get_fits("Epoch_6")
#unbin['usum'].name = 'flx'
#spec_bin binning - also calculates and includes Ptot and Perr 
specbin = Bin_data(unbin, old_binSize= 4, binSize = 20, z = 0.00483)

#binning as descibbed in slack
#oldbin = old_Bin_data(unbin, old_binSize = 4, binSize = 20, z = 0.00483)

fig, (ax0, ax1, ax2, ax3, ax4) = plt.subplots(5, sharex = True, figsize=(20,14))
fig.suptitle("Binning Comparison - 2010jl", fontsize=24, y=.95)
ax0.plot(specbin['wave'], specbin['flx'], c = 'b', label = 'Spectra_binning'); ax0.set_ylabel('Flux') 
#ax0.step(oldbin['wave'], oldbin['flx'], c = 'r', label = 'Hand_calculated_binning')
ax1.step(specbin['wave'], specbin['qsum'], c = 'b'); ax1.set_ylabel('Qsum')  
#ax1.step(oldbin['wave'], oldbin['qsum'], c = 'r')
ax2.step(specbin['wave'], specbin['p'], c = 'b') ; ax2.set_ylabel('%P') 
#ax2.plot(oldbin['wave'], oldbin['p'], c = 'r')
ax3.step(specbin['wave'], specbin['q'], c = 'b'); ax3.set_ylabel('%Q')  
#ax3.plot(oldbin['wave'], oldbin['qerr'], c = 'r')
ax4.step(specbin['wave'], specbin['u'], c = 'b') ; ax4.set_ylabel('%U') 
#ax4.plot(oldbin['wave'], oldbin['uerr'], c = 'r')
fig.legend()
"""

#---------------Get data for UGs------------
"""
#for 2014ad
unbin3 = get_txt_data("SN_2014ad/epoch3/", "sn2014ad_2014-05-22comb.", "2014ad_2014-05-22all_data")
#unbin3 = get_txt_data("SN_2014ad/epoch2/", "sn2014ad_2014-04-26comb.", "2014ad_2014-04-26all")
#unbin3 = get_txt_data("SN_2014ad/epoch3/", "sn2014ad_2014-05-22comb.", "2014ad_2014-05-22all")
epoch3 = Bin_data(unbin3, 4, 30, z = 0.005723, new_filename = "SN_2014ad/epoch3/2014ad_2014-05-22all_binned30" )

print(epoch3['q'])
fig, (ax0, ax1, ax2) = plt.subplots(3, sharex = True, figsize=(10,8))
ax2.step(epoch3['wave'], epoch3['p'], c = 'k')
ax0.step(epoch3['wave'], epoch3['q'], c ='k')
ax1.step(epoch3['wave'], epoch3['u'], c = 'k')
ax0.set_ylabel('Q');ax0.plot(unbin3['wave'], unbin3['q']*100); ax0.set_ylim(-2, 2)
ax1.set_ylabel('U');ax1.plot(unbin3['wave'], unbin3['u']*100); ax1.set_ylim(-2, 2)
ax2.set_ylabel('P');ax2.plot(unbin3['wave'], unbin3['p']*100); ax2.set_ylim(-2, 2)
"""
"""
#for 2012au - did not redo and upload on onedrive yet
unbin = get_txt_data("Epoch_1/all_comb/", "sn2012au_allcomb.", "SN_2012au-2012-03-21all_data")
binned = Bin_data(unbin, 4, 30, new_filename="Epoch_1/all_comb/SN_2012au-2012-03-21all_binned30")
"""
"""
#for M12045
unbin1 = get_txt_data("SN_M12045/epoch7/", "M12045_2015-04-27comb.", "M12045_2015-04-27all_data")
epoch1_20 = Bin_data(unbin1, 2.5, 10, z = 0.001891, new_filename="SN_M12045/epoch7/M12045_2015-04-27all_binned10")

fig, (ax0, ax1, ax2, ax3, ax4) = plt.subplots(5, sharex = True, figsize=(20,14))
fig.suptitle("Epoch 7 M12045", fontsize=24, y=.95)
ax0.plot(epoch1_20['wave'], epoch1_20['flx']*10**15, c = 'b', label = 'Spectra_binning'); ax0.set_ylabel('Flux') 
ax0.plot(epoch1_20['wave'], epoch1_20['qsum']/10**4.7, c = 'r'); 
ax1.step(epoch1_20['wave'], epoch1_20['p'], c = 'b') ; ax1.set_ylabel('%P') 
ax2.step(epoch1_20['wave'], epoch1_20['q'], c = 'b'); ax2.set_ylabel('%Q')  
ax3.step(epoch1_20['wave'], epoch1_20['u'], c = 'b') ; ax3.set_ylabel('%U') 
ax4.step(epoch1_20['wave'], PA(epoch1_20['p'], epoch1_20['u']), c = 'b') ; ax4.set_ylabel('PA') 
fig.legend()
"""

#for ASASSN-16fp - last binned 2/9/24, binned epoch 6 on 5/29/24, still need epoch 1 (other MMT epoch)
"""
unbin = get_txt_data("SN_ASASSN16fp/epoch6/", "ASASSN16fp_2016-11-06comb.","ASASSN16fp_2016-06-26all" )
binned = Bin_data(unbin, 4, 30, z = 0.00365, new_filename= "SN_ASASSN16fp/epoch6/ASASSN16fp_2016-11-06all_binned30")
"""
#--------------Bin Data for Manisha SN 2023 ixf----------
"""
#looking at fits data: 
#hdulist = pyfits.open("SN2023ixf_polcombined/jun18/SN2023ixf_allcomb_23jun18.qsig.fits")
#print(hdulist[0].header)
#tbdata = hdulist[0].data
#print(tbdata)    

binsize = 20
epochs = ['may31', 'jun01', 'jun02', 'jun03', 'jun04', 'jun05']
for folder in epochs:
    data = get_fits("SN2023ixf_polcombined/"+folder, z = 0.000804) #this is where deredshift should happen (it was not here when I binned for Manisha so hers is not deredshifted)
    data.rename_column('usum', 'flx')
    binned = Bin_data(data, 4, binsize)
    binned = get_pa_column(binned)
    new_filename = "SN2023ixf_polcombined/SN2023_binned/SN2023ixf_"+folder+"_b20alldat.txt"
    binned.write(new_filename+".txt", format = 'ascii')
"""