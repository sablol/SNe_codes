# -*- coding: utf-8 -*-
"""
Created on Fri Oct 15 11:47:14 2021

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

#open, combine and bin data
data1 = get_txtFITS("Epoch_1/all_comb", "sn2012au_allcomb.qsum.txt", "sn2012au_allcomb.q.txt", "sn2012au_allcomb.qsig.txt", "sn2012au_allcomb.flx.txt", "sn2012au_allcomb.u.txt", "sn2012au_allcomb.usig.txt")
epoch1_20 = Bin_data(data1, 4, 20)
data2 = get_txtFITS("Epoch_2", "sn2012au_allcomb.qsum.txt", "sn2012au_allcomb.q.txt", "sn2012au_allcomb.qsig.txt", "sn2012au_allcomb.flx.txt", "sn2012au_allcomb.u.txt", "sn2012au_allcomb.usig.txt")
epoch2_20 = Bin_data(data2, 2.5, 20)
data3 = get_txtFITS("Epoch_3", "sn2012au_allcomb.qsum.txt", "sn2012au_allcomb.q.txt", "sn2012au_allcomb.qsig.txt", "sn2012au_allcomb.flx.txt", "sn2012au_allcomb.u.txt", "sn2012au_allcomb.usig.txt")
epoch3_20 = Bin_data(data3, 4, 20)
data4 = get_txtFITS("Epoch_4", "sn2012au_allcomb.qsum.txt", "sn2012au_allcomb.q.txt", "sn2012au_allcomb.qsig.txt", "sn2012au_allcomb.flx.FIXED.txt", "sn2012au_allcomb.u.txt", "sn2012au_allcomb.usig.txt")
epoch4_20 = Bin_data(data4, 4, 20); epoch4_20 = epoch4_20[10:]
data5 = get_txtFITS("Epoch_5", "sn2012au_allcomb.qsum.txt", "sn2012au_allcomb.q.txt", "sn2012au_allcomb.qsig.txt", "sn2012au_allcomb.flx.txt", "sn2012au_allcomb.u.txt", "sn2012au_allcomb.usig.txt")
epoch5_20 = Bin_data(data5, 4, 20)
data6 = get_fits("Epoch_6")
data6['usum'].name = 'flx'
epoch6_20 = Bin_data(data6, 2.5, 20)


e1 = plot_line_velspace(epoch1_20, 5876, 20000, 10**-5, "Epoch1(Day0-7)")
#e2 = plot_line_velspace(epoch2_20, 5876, 20000, 10**-5.7, "Epoch2(Day26)")
#e3 = plot_line_velspace(epoch3_20, 5876, 20000, 10**-4.5, "Epoch3(Day35-40)")
#e4 = plot_line_velspace(epoch4_20, 5876, 20000, 10**-4.3, "Epoch4(Day57-67)")
#e5 = plot_line_velspace(epoch5_20, 5876, 20000, 10**-4, "Epoch5(Day85-90)")
#e6 = plot_line_velspace(epoch6_20, 5876, 20000, 10**-4.7, "Epoch6(Day295)")



min_wave = np.round(get_min(e1, 'p')[0])
print(min_wave)
print(get_velocity(5876, min_wave ))
        
He = [5876, min_wave ]; He_lab = ["He I (5876)"]; He_c = ["b", 'm']
plot_p_single(He, He_lab, He_c, data1, epoch1_20, 20, 10**-4.6, "Epoch_1")

#make_gif(["He_5876/Epoch_Polarization_Mins/1.png", "He_5876/Epoch_Polarization_Mins/2.png", "He_5876/Epoch_Polarization_Mins/3.png", "He_5876/Epoch_Polarization_Mins/4.png", "He_5876/Epoch_Polarization_Mins/5.png", "He_5876/Epoch_Polarization_Mins/6.png"], "He_5876_pol_min.GIF", 2)


