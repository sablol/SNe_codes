# -*- coding: utf-8 -*-
"""
Created on Wed Jan 20 13:04:03 2021

@author: Rachel
"""
####Edited this from Rachels version on one drive refer back to that if there are issues
#this one is designed to only worry about the flux data


#%% Import Librries
from astropy.io import ascii,fits
from astropy.table import Table,hstack,vstack,Column
import numpy as np
import os
import matplotlib.pyplot as plt
import glob
#from lines_pol_module import LinePol
import scipy.stats as stats
from scipy import interpolate
import sympy as sym
from sympy.solvers import solve

#%% Read in the spectra and phase data
dir = 'D://Documents//DU//WR071'
# phases = ascii.read(dir + '//' + 'WR71_Phases.txt', data_start=1)
 
#-----------TO DO-----------------------------------
#set blue continuum and red continuum interval values 
#when choosing 3 points make sure to get mins not maxes for absorption
#figure out where variance value comes from (data or fit?)

#-------------get SN flux data-------------------------------------------------
import matplotlib.pyplot as plt #this is a package that helps to plot things 
from astropy.io import ascii #this is a package to help open the file

data = ascii.read("SN2014ad_2014-05-22comb.flx.txt") #reading in the data as lon as the datafile is in the same folder as this script file

wavelengths = [] #empty list to fill in with wavelength values
fluxs = [] #empty list to fill in with flux values

for i in range(len(data)): #loops through each line in file and 
    wavelengths.append(data[i][0]) #appends wavelength value in that in to wavelength list
    fluxs.append(data[i][1]) #appends flux value in that line to flux list 

#plot wavlength list vs. flux list
plt.figure(figsize=(10, 8)) #creates a figure and specifies the size (base, hight)
plt.plot(wavelengths, fluxs, c = 'k', linewidth = 1)#plots data given (x axis data, y axis data, c = color, linewidth = thickness of line)
plt.xlabel('wavelength(Å)')#label x axis 
plt.ylabel('relative flux')#label y axis
plt.title('SN2018_gv-epoch5')#add title 


#Identify regions --- Sodium MW doublet
bluewavbeg = 5870
bluewavfin = 5885
linestart = 5887
linend = 5898
redwavbeg = 5900
redwavfin = 5910
"""
#for i in data_list:
# Identify line region and sort       
#    maximum = np.amax(stack_table['Stokes I'][lower_wave:upper_wave])
linereg = []

for i in range(len(array[:])):
    if linestart <= array[i][0] and  linend >= array[i][0]:
        linereg.append(array[i])

linereg = np.array(linereg)

fullreg = []

for i in range(len(array[:])):
    if bluewavbeg  <= array[i][0] and redwavfin >= array[i][0]:
        fullreg.append(array[i]) 

fullreg = np.array(fullreg)

iSorted = linereg[linereg[:,1].argsort()[::-1]]

# Interpolate Continuum
#averages blue and red side by drawing a line between mean value on either side
bluereg = np.arange(bluewavbeg, bluewavfin, 1) #wave range of continuum on blue side of feature 
redreg = np.arange(redwavbeg,redwavfin,1)    
inter = np.hstack((bluereg, redreg)) 
contreg = []

for i in range(len(array[:])):
    if bluewavbeg <= array[i][0] and bluewavfin >= array[i][0]:
        contreg.append(array[i])

for i in range(len(array[:])):
    if redwavbeg <= array[i][0] and redwavfin >= array[i][0]:
        contreg.append(array[i])

contreg = np.array(contreg)
 
iglobalvar.append((sum(contreg[:,2])/len(contreg[:,2]))) #this if the flux
#qglobalvar.append(sum(contreg[:,4])/len(contreg[:,4]))
#uglobalvar.append(sum(contreg[:,6])/len(contreg[:,6]))
#pglobalvar.append(sum(contreg[:,8])/len(contreg[:,8]))
        
ifit = interpolate.interp1d(contreg[:,0], contreg[:,1], kind='linear', axis= -1, copy=True, bounds_error=None, assume_sorted=False, fill_value='extrapolate')
#qfit = interpolate.interp1d(contreg[:,0], contreg[:,3], kind='linear', axis= -1, copy=True, bounds_error=None, assume_sorted=False, fill_value='extrapolate')
#ufit = interpolate.interp1d(contreg[:,0], contreg[:,5], kind='linear', axis= -1, copy=True, bounds_error=None, assume_sorted=False, fill_value='extrapolate')
#pfit = interpolate.interp1d(contreg[:,0], contreg[:,7], kind='linear', axis= -1, copy=True, bounds_error=None, assume_sorted=False, fill_value='extrapolate')

# Preparing to subtract continuum from spectra
ifunc = ifit(fullreg[:,0])
#qfunc = qfit(fullreg[:,0])
#ufunc = ufit(fullreg[:,0])
#pfunc = pfit(fullreg[:,0])


ifunc = np.vstack((ifunc,fullreg[:,0],fullreg[:,2]))
#qfunc = np.vstack((qfunc,fullreg[:,0], fullreg[:,4]))
#ufunc = np.vstack((ufunc,fullreg[:,0], fullreg[:,6]))
#pfunc = np.vstack((pfunc,fullreg[:,0], fullreg[:,8]))


ifunc = np.transpose(ifunc)
#qfunc = np.transpose(qfunc)
#ufunc = np.transpose(ufunc)
#pfunc = np.transpose(pfunc)


icontsubspec = []
 
# plt.plot(ifunc[:,1], ifunc[:,0])
# plt.plot(fullreg[:,0], fullreg[:,1])
# plt.show()

iroi = [] 

for i in range(len(ifunc)):
    if linestart <= ifunc[i][1] and  linend >= ifunc[i][1]:
        iroi.append((ifunc[i,0], ifunc[i,1], ifunc[i,2]))

iroi = np.array(iroi)

        
# Continuum Subtraction

for i in range(len(iroi[:])):
    icontsubspec.append(abs(linereg[i,1] - iroi[i]))

icontsubspec = np.array(icontsubspec)
  
#Three largest points 
iSorted = icontsubspec[np.argsort(icontsubspec[:,0])[::-1]] #icontsubspec if flux data and found peak points of line comparing each point to the one next to it
    
# Find max stokes values and errors within line region 
iwave = [iSorted[0][1], iSorted[1][1], iSorted[2][1]]
ival = [iSorted[0][0], iSorted[1][0], iSorted[2][0]]   
############NOTE!!!!!!!!! Variances are from fits files with absolutely no propogation!!!!!!!!!!!!##################  
ivar = [iSorted[0][2], iSorted[1][2], iSorted[2][2]]
   
ivaravg = sum(ivar)/len(ivar)

iweight = []
 
for i in ivar:
    ivar = abs(i)
    iweight.append(1/ivar)
    
# Quadratic fit to weighted top 3 flux values
imodel = np.polyfit(iwave,ival, 2, w = iweight) #2 mean quadratic fit (1 = line, 3 = cubical, ect)
if iSorted[2][1] > iSorted[0][1]:
    iSorted[2][1] = iSorted[2][1]+1
    iSorted[0][1] = iSorted[0][1]-1
else:
    iSorted[2][1] = iSorted[2][1]-1
    iSorted[0][1] = iSorted[0][1]+1
    
# Find stokes value of peak
ipolyline = np.linspace(iSorted[2][0], iSorted[0][0], 500)
   
ieq = (ipolyline**2)*imodel[0] + (ipolyline)*imodel[1] + imodel[2] #peak value of flux in fit
   
idiff = -imodel[1]/(2*imodel[0]) #solving quadratic for x value where function is at its peak (should be local min for me)

#    plt.scatter(wave, ival)
#    plt.plot(polyline, eq)    
#    plt.figure()

#------------------------------------------------------------------ 
#------------------------------------- I data
wave_in =(linereg[:,0])
I_flux_in = (linereg[:,1])
I_flux_var =(linereg[:,2]) # I don't have this column in my data, not sure where this comes from 
I_flux_err = np.sqrt(linereg[:,2])

#-------------------------------------- Integration method (Trapezoid rule), Chris 
#-------------------------------------- let i, q, u be the long sides of the trapezoids in units of flux/Angstroms
#-------------------------------------- width is the wavelength width, Angstroms (in this case, all of the wavelengths are 
#---------------------------------------1 A apart)
#-------------------------------------- calculating the total Intensity within each 1 A bin

I_total_list = []
width = 4 #angstrom resolution 
for i1, i2 in zip(I_flux_in[0::], I_flux_in[1::]): #i1 and i2 are heights of each side of trapazoid
    I_total_list.append(width*(i1 + i2)/2) #calculate are in trapezoid and append to list 

#----------------------------------------Total I    

totI=np.cumsum(I_total_list) #sum up all areas found in trapezoids [a1, a1+a2, a1+a2+a3...]
IINT = np.max(totI) #finds actual sum of areas (largest value in list above)
#print("Summing the trapezoids area up to find the total area yeilds:")
#print()
#print(totalarea)
#print()

# Calculate width of line region in I, Q, and U

iwidth = abs(IINT)/np.max(ieq) #calculate width of line that trapezoid fit into
lineWidthI.append(iwidth)

 #not sure about what to do with var stuff   
iglobalvar = np.array(iglobalvar)

globalvariances = Table([iglobalvar], names = ['Ivar'])
"""

