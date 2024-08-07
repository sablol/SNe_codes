# -*- coding: utf-8 -*-
"""
Created on Mon Feb 26 09:23:15 2024

@author: sabri
"""

import numpy as np
import matplotlib.pyplot as plt
import math as math
import random
from astropy.table import Table, setdiff
import time


def make_grid(nphot = [1003, 1003], lim = [6,6]):
    #Reproduce Maund/Reilly Model
    #start = time.time() #record start time 
    #create grid of photons (identify each photon(pts) by x,y coordinate)
    nx, ny = nphot[0], nphot[1] #1003, 1003 #number of photons(pts) in x and y directions, for a center pt nx = whole Number*xlim +1 (same for ny)
    xlim, ylim = lim[0], lim[1] #grid size limit on either side of center (arbirary size (using units of AU) to fit photosphere)
    x = np.linspace(-xlim, xlim, nx) #values that grid will cover (using units of AU)
    y = np.linspace(-ylim, ylim, ny)
    midx = int((len(x) - 1)/2); midy = int((len(y) - 1)/2)
    centerx = x[midx]; centery = y[midy] #specify center
    #end = time.time(); time_sec = (end-start); print("Func: make_grid run time:", (time_sec), "seconds") 
    return(x, y, centerx, centery)

def make_photosphere(grid_pts, center, phot_axis, phot_rotation, contour = 'off'):
    start = time.time() #record start time
    x, y = grid_pts
    ella, ellb = phot_axis #photosphere major axis smaller than bounds of grid
    ell_thetas = np.linspace(0, 2*math.pi, 100) #angles - theta = 0 = +x , theta = pi = -x, theta = pi/2 = +y, theta = 3pi/2 = -y
    rotation = 2*math.pi-np.radians(phot_rotation) #photosphere major axis rotated 119degrees from Reilly 2016 rotates from x = 0 degrees (360-119 degrees but runs in radian)
    Ell = np.array([ella*np.cos(ell_thetas) , ellb*np.sin(ell_thetas)]) #ellipse radius pts in [x, y] array
    R_rot = np.array([[np.cos(rotation) , -np.sin(rotation)], [np.sin(rotation) , np.cos(rotation)]]) #2-D rotation matrix
    Ell_rot = np.zeros((2,Ell.shape[1])) #empty [x,y] array for rotated radius pts 
    for i in range(Ell.shape[1]):
        Ell_rot[:,i] = np.dot(R_rot,Ell[:,i]) #rotate ellipse by dot product of rotation matrix and original ellipse x,y coords array
    
    centerx, centery = center
    phot_xcoords = centerx+Ell_rot[0,:] #photosphere ellispe coordinates 
    phot_ycoords = centery+Ell_rot[1,:]
    phot_coords = (phot_xcoords, phot_ycoords)
    
    #hold p, PA, I and x,y coords for each photon (x, y pt) in photosphere
    p_in_phot = []
    I_in_phot = []
    x_in_phot = []
    y_in_phot = []
    PA_in_phot = []
    zp = []; zI =[] #arrays for %p and I values for each pt in grid (z dimention for countour plot)
    #assign p values to pts in grid based on radius
    for rx in x:
        dummy = [] #dummy temp lists for countour plots 
        dum = []
        for ry in y: 
            r = np.sqrt(rx**2 + ry**2) #calculate radius at current pt
            theta = math.atan2(ry, rx) #calculate theta at current pt
            r_phot = np.sqrt((centerx+ella*np.cos(theta-rotation))**2 + (centery+ellb*np.sin(theta-rotation))**2) #calculate radius at pt on photosphere at same angle as pt in grid (account for rotated photosphere) 
            if r < r_phot:
                rphotx = r_phot*np.cos(theta); rphoty = r_phot*np.sin(theta) #calculate corresponding x,y coord on photosphere radius to rx,ry coord for pt (photon) in photosphere (r_phot already accounts for rotated photosphere)
                p = 0.15*(r/r_phot) #calculate %p equation from Maund 2010 (constant for %15 of light at limb polarized)
                I = 0.5*(1-np.sqrt(1-((r/r_phot)**2))) #I = Ir/I_0, the ratio of the lumonosity, using k = 0.5 for limb darkening from Maund 2010
                p_in_phot.append(p) #keep track of %p for pts in the photosphere
                I_in_phot.append(I) #keep track of Ir/I_0 for pts in the photosphere 
                x_in_phot.append(rx); y_in_phot.append(ry) #keep track of x, y coords. for pts in photosphere
                #calculate PA's used for model (16.5% parallel to photosphere, the rest randomly assigned)
                grn = random.uniform(0, 1) #generate random number to identify each pt (or photon) in phot by
        
                if grn > 0.165: # when random number > 0.165 assign random angle (between 0 and 180) to list
                    #grn2 = random.uniform(0, 1) #generate new random number to calculate random PA 
                    #PA_in_phot.append(180*grn2) # calculate and assign random angle (between 0 and 180) to list 
                    PA_in_phot.append(random.uniform(0, 180))#generate random angle between 0 and 180
                else:
                    #if random number < 0.165 append PA parallel to limb - this should cover that 16.5% of PA's be parallel to the limb Reilly 2016
                    #find slope for pts in ellipse according to https://www.anirdesh.com/math/algebra/ellipse-tangents.php
                    m = (ellb*rphotx)/(ella*np.sqrt(ella**2 - rphotx**2)) #calculate slope at corresponding location on ellipse 
                    if (rx <= 0 and ry >= 0) or (rx >= 0 and ry <= 0):
                        m = np.abs(m) #pts with positive slope (for ellipse rotated clockwise from major axis at x=0):
                    else: 
                        #negative slope
                        if m > 0: 
                            m = m*-1
                    angle = math.degrees(math.atan(m))#get angle of slope
                    if angle < 0:
                        angle = angle + 180 #add 180 to negative so in QU plane (0-180)
                    PA_in_phot.append(angle) #assigns every pt (photon) a PA that is the angle of slope in degrees
            if contour == 'contour':
                ###calc p here for contour plot otherwise saves time to do it only for pts in phot
                p = 0.15*(r/r_phot); #calc p here for countour plot
                I = 0.5*(1-np.sqrt(1-((r/r_phot)**2))); #calc I here for countour plot
                dummy.append(p) #(for countour plot)
                dum.append(I)
        if contour == 'contour':
            zp.append(dummy) #(for countour plot)
            zI.append(dum)
    #calculate q*I and u*I for pts in photosphere (do all math in intrinsic units)
    qI = (p_in_phot*np.cos(2*np.radians(PA_in_phot)))*I_in_phot  
    uI = (p_in_phot*np.sin(2*np.radians(PA_in_phot)))*I_in_phot
    
    #Create Table of values for pts (photons) in photosphere
    In_Phot = Table([x_in_phot, y_in_phot, p_in_phot, PA_in_phot, qI, uI, I_in_phot], names = ['x', 'y', 'p', 'pa', 'qI', 'uI', 'I'])
    
    end = time.time(); time_sec = (end-start); print("Func: make_photosphere run time:", (time_sec), "seconds")  
    return(phot_coords, In_Phot, zp, zI)


def make_clump(center, clump_axis, clump_angle, In_Phot_dat):
    start = time.time() #record start time
    #ella_clump, ellb_clump = clump_ella, clump_ellb #3.5, 2.485 #ellipse radius ella = x, ellb = y so clump_ella = .88*phot_ella and ellb/ella = 0.71 from Reilly 2016
    ella_clump, ellb_clump = clump_axis
    ell_thetas_clump = np.linspace(0, 2*math.pi, 100) #angles - theta = 0 = +x , theta = pi = -x, theta = pi/2 = +y, theta = 3pi/2 = -y
    rotation_clump = 2*math.pi-np.radians(clump_angle) #clump major axis rotated 160 degrees from Reilly 2016 rotates from x = 0 degrees 
    Ell_clump = np.array([ella_clump*np.cos(ell_thetas_clump) , ellb_clump*np.sin(ell_thetas_clump)]) #clump radius x,y coords array 
    R_rot_clump = np.array([[np.cos(rotation_clump) , -np.sin(rotation_clump)], [np.sin(rotation_clump) , np.cos(rotation_clump)]]) #2-D rotation matrix
    Ell_rot_clump = np.zeros((2, Ell_clump.shape[1])) #empty rotated clump radius x,y coords array 
    for i in range(Ell_clump.shape[1]):
        Ell_rot_clump[:,i] = np.dot(R_rot_clump, Ell_clump[:,i]) #rotate x,y coords and fill in empty array
    centerx, centery = center
    clump_xcoords = centerx+Ell_rot_clump[0,:] #rotated clump radius x coords, displaced from center 
    clump_ycoords = centery+Ell_rot_clump[1,:] #rotated clump radius x coords, displaced from center
    clump_coords = (clump_xcoords, clump_ycoords)
    
    #Isolate data for pts in clump going through only photosphere pts (ie clump cannot extend past the bounds of the photosphere)
    In_Clump = Table(names = ['x', 'y', 'p', 'pa', 'qI', 'uI', 'I']) #creat table for pts within clump radius
    for ph_dat in In_Phot_dat:
        r = np.sqrt(ph_dat['x']**2 + ph_dat['y']**2) #calculate radius at current pt in photosphere
        theta_c = math.atan2(ph_dat['y'], ph_dat['x']) - rotation_clump #calculate theta at current pt (account for clump rotation)
        r_clump = np.sqrt((centerx+ella_clump*np.cos(theta_c))**2 + (centery+ellb_clump*np.sin(theta_c))**2) #calculate radius of clump at same angle as current pt in photosphere
        if r < r_clump:
           In_Clump.add_row(ph_dat) #add values for pts inside clump radius into table
    
    end = time.time(); time_sec = (end-start); print("Func: make_clump run time:", (time_sec), "seconds")       
    return(clump_coords, In_Clump)


def get_InPhot_NoClump(In_Phot_dat, In_Clump_dat):
    start = time.time() #record start time
    #Isolate data in photosphere but not obscured by clump ------------------------------------------------------------------------
    Phot_noClump = setdiff(In_Phot_dat, In_Clump_dat)
    #calculate %p, PA in photosphere region not obscured by clump 
    #*NOTE* this value changes because the PA angle of each pt in grid is randomly assigned (except for the 16.5% that is parallel to the photosphere)
    Tot_P = (np.sqrt(np.sum(Phot_noClump['qI'])**2 + np.sum(Phot_noClump['uI'])**2))/np.sum(Phot_noClump['I']) #sum total q and u values by vector decomp. & divide by total I to get total p  
    qsum = np.sum(Phot_noClump['qI'])/np.sum(Phot_noClump['I']) #get q and u values to calculate new PAs
    usum = np.sum(Phot_noClump['uI'])/np.sum(Phot_noClump['I'])
    Tot_PA = np.degrees(.5*np.arctan2(usum, qsum)) #calculate PA by q and u vals 
    if Tot_PA < 0:
        Tot_PA = Tot_PA + 180
    
    end = time.time(); time_sec = (end-start); print("Func: get_InPhot_NoClump run time:", (time_sec), "seconds")  
    return(Tot_P*100, Tot_PA)

def MR_2Dmodel(nphot:list, photsph_ratio:float, photsph_angle:float, clumpphot_frac:float, clump_ratio:float, clump_angle:float, plot = 'off'):
    x, y, centx, centy = make_grid(nphot)  
    
    ella_phot = np.max(x)- 1 #photosphere major axis smaller than bounds of grid
    ellb_phot = ella_phot*photsph_ratio
    photsph_coords, InPhot_vals, zp, zI = make_photosphere(grid_pts = (x, y), phot_axis = (ella_phot, ellb_phot), phot_rotation = photsph_angle, center = [centx, centy], contour = plot)
    
    ella_clump = ella_phot*clumpphot_frac
    ellb_clump = ella_clump*clump_ratio
    clump_coords, InClump_vals = make_clump(center = [centx, centy], clump_axis = (ella_clump, ellb_clump), clump_angle = clump_angle, In_Phot_dat=InPhot_vals)
    
    Tot_P, Tot_PA = get_InPhot_NoClump(In_Phot_dat = InPhot_vals, In_Clump_dat = InClump_vals)
    
    if plot == 'on':
        ####Plot model set up to check (only for few number pf photons (nx, ny)) 
        #make meshgrid from photon pts 
        xp, yp = np.meshgrid(x, y) #x, y values for individual grid pts
        fig, ax = plt.subplots(1)
        ax.set_xticks(np.arange(np.min(x), np.max(x)+1, 2)); ax.set_yticks(np.arange(np.min(y), np.max(y)+1, 2)); ax.axis('square')
        ax.plot(xp, yp, marker='o', color='k', linestyle='none') #plot grid pts (individual photons)
        ax.plot(centx, centy, marker='o', color='r', linestyle='none') #plot center
        ax.plot(InPhot_vals['x'], InPhot_vals['y'], marker='o', color='m', linestyle='none') #plot pts (photons) in photosphere 
        ax.plot(InClump_vals['x'], InClump_vals['y'], marker='o', color='lime', linestyle='none') #plot pts (photons) in clump 
        #ax.plot(Phot_noClump['x'], Phot_noClump['y'], marker='o', color='lightblue', linestyle='none') #plot pts (photons) in photosphere but not in clump 
        ax.plot(photsph_coords[0], photsph_coords[1], color = 'r') #plot photosphere radius
        ax.plot(clump_coords[0] , clump_coords[1], 'y')#plot clump radius
      
    if plot == 'contour':
        ###countour plot of I or P values by radius 
        xp, yp = np.meshgrid(x, y) #x, y values for individual grid pts
        xs, ys = np.meshgrid(x, y, sparse = True)#x, y values for countour plot
        #plot %p and I values as countour plot
        fig, ax = plt.subplots(1, 2)
        ax[0].set_xticks(np.arange(np.min(x), np.max(x)+1, 2)); ax[0].set_yticks(np.arange(np.min(y), np.max(y)+1, 2)); ax[0].axis('square')
        ax[0].plot(xp, yp, marker='o', color='k', linestyle='none') #plot grid pts (individual photons)
        ax[0].plot(centx, centy, marker='o', color='r', linestyle='none') #plot center
        contP = ax[0].contourf(x, y, zp)
        cbar = fig.colorbar(contP, orientation='horizontal'); cbar.set_label('% P')#cbar = fig.colorbar(plot, orientation='horizontal'); cbar.set_label('I')

        ax[1].set_xticks(np.arange(np.min(x), np.max(x)+1, 2)); ax[1].set_yticks(np.arange(np.min(y), np.max(y)+1, 2)); ax[1].axis('square')
        ax[1].plot(xp, yp, marker='o', color='k', linestyle='none') #plot grid pts (individual photons)
        ax[1].plot(centx, centy, marker='o', color='r', linestyle='none') #plot center
        contI = ax[1].contourf(x, y, zI)
        cbar = fig.colorbar(contI, orientation='horizontal'); cbar.set_label('I')

    print('Model ran with ', nphot[0]*nphot[1], 'total photons')
    print('Photosphere Ellipticity: ', ellb_phot/ella_phot, ', rotated by: ', photsph_angle, '°')
    print('Clump coverage: ', clumpphot_frac*100, '% of Photosphere')
    print('Clump Ellipticity: ', ellb_clump/ella_clump, ', rotated by: ', clump_angle, '°')
    print('Returns: Total P(%): ', np.round(Tot_P, 2), ' Total PA(°): ', np.around(Tot_PA, 2), ' Total Photons: ', nphot[0]*nphot[1])    

    return(Tot_P, Tot_PA, nphot[0]*nphot[1])


#testing multiple model runs with same inputs (~1/2million photons) to plot results and see if the P and PA stabilizes 
t = Table(names = ['P(%)', 'PA(°)'])
Ptots = []; PAtots = []; Nphotons = []
for i in np.arange(0, 25):
    #(721*721 = 519841 total photons) (nx*ny = total photosn but nx & ny must be in increments of x*6+1)
    ixf_mod = MR_2Dmodel(nphot = [901, 901], photsph_ratio = 0.83, photsph_angle = 165, clumpphot_frac = 0.88, clump_ratio = 0.71, clump_angle = 265, plot = 'off')
    Ptots.append(ixf_mod[0]); PAtots.append(ixf_mod[1]); Nphotons.append(ixf_mod[2])
    t.add_row((ixf_mod[0], ixf_mod[1]))

t.write("MR_nphot811801.txt", format = 'ascii')    
fig, ax = plt.subplots(1,1, sharex = True);plt.subplots_adjust(hspace = 0)
ax.set_xlabel('PA (°)'); ax.set_ylabel('Total P(%)')
ax.plot(PAtots, Ptots, marker = 'o', c = 'grey', mfc = 'm')


"""
iter = 1
t = Table(names = ['P(%)', 'PA(°)'])
while iter <= 10:
    start = time.time() #record start time
    out = MR_2Dmodel(nphot = [1003, 1003], photsph_ratio = 0.91, photsph_angle = 119, clumpphot_frac = 0.88, clump_ratio = 0.71, clump_angle = 160, plot = 'off')
    end = time.time(); time_min = (end-start)/60; print("Func: get_InPhot_NoClump run time:", (time_min), "minutes")  
    t.add_row((out[0], out[1]))
    print('on round', iter)
    iter += 1
t.write("Maund_Reilly_Totals.txt", format = 'ascii')
"""
"""
t = Table(names = ['P(%)', 'PA(°)'])
out = MR_2Dmodel(nphot = [1003, 1003], photsph_ratio = 0.91, photsph_angle = 119, clumpphot_frac = 0.88, clump_ratio = 0.71, clump_angle = 160, plot = 'on')
t.add_row((out[0], out[1]))
t.write("Maund_Reilly_Totals.txt", format = 'ascii')
"""
"""
def MR_2Dmodel(phot_ella, phot_ellb, phot_rotation, clump_ella, clump_ellb, clump_rotation, plot = 'off'):
        
    #For photosphere ---------------------------------------------------------------------------------------------------------
    #plot ellipse photosphere radius
    
    
    #For obscuring clump ----------------------------------------------------------------------------------------------------
        
    #Isolate data in photosphere but not obscured by clump ------------------------------------------------------------------------
    Phot_noClump = setdiff(In_Phot, In_Clump)
    #calculate %p, PA in photosphere region not obscured by clump 
    #*NOTE* this value changes because the PA angle of each pt in grid is randomly assigned (except for the 16.5% that is parallel to the photosphere)
    Tot_P = (np.sqrt(np.sum(Phot_noClump['qI'])**2 + np.sum(Phot_noClump['uI'])**2))/np.sum(Phot_noClump['I']) #sum total q and u values by vector decomp. & divide by total I to get total p  
    qsum = np.sum(Phot_noClump['qI'])/np.sum(Phot_noClump['I']) #get q and u values to calculate new PAs
    usum = np.sum(Phot_noClump['uI'])/np.sum(Phot_noClump['I'])
    Tot_PA = np.degrees(.5*np.arctan2(usum, qsum)) #calculate PA by q and u vals 
    if Tot_PA < 0:
        Tot_PA = Tot_PA + 180
    #print("Total % p = " + str(Tot_P*100) + "Position Angle =" + str(Tot_PA)) #Values from Reilly are p = 0.7 +/- 0.5, theta = 70 +/- 3
    
    if plot == 'on':
        plot_mod()
        
    
    end = time.time(); time_min = (end-start)/60; print("Model run time:", (time_min), "minutes")
    print('')
    
    return(Tot_P*100, Tot_PA, nx*ny)
"""
#print(MR_2Dmodel(1008, 1008, 5, 4.55, 119, 3.5, 2.485, 160))
#print(MR_2Dmodel(1008, 1008, 5, 4.55, 241, 3.5, 2.485, 200)) #check other angle rotation
#print(MR_2Dmodel(1008, 1008, 4, 2, 119, 3, 2, 160)) #check other ellipse ratios
"""
#testing up to 1million photons to plot results and see at what nphoton the P and PA stabilizes 
Ptots = []; PAtots = []; Nphotons = []
for i in np.arange(36, 1008, 36):
    #(1008*1008) (nx*ny = total photosn but nx & ny must be in increments of x*6+1)
    nx = ny = i+1
    #out = MR_2Dmodel(nx, ny, 5, 4.55, 119, 3.5, 2.485, 160)
    Ptots.append(out[0]); PAtots.append(out[1]); Nphotons.append(out[2])
    
fig, ax = plt.subplots(2,1, sharex = True);plt.subplots_adjust(hspace = 0)
ax[1].set_xlabel('Number of Photons Ran')
ax[0].plot(Nphotons, Ptots, marker = 'o', c = 'grey', mfc = 'm'); ax[0].set_ylabel('Total P(%)')
ax[1].plot(Nphotons, PAtots, marker = 'o', c ='grey', mfc = 'c'); ax[1].set_ylabel('Total PA(°)')
""" 
"""  
def MR_2Dmodel(nx, ny, phot_ella, phot_ellb, phot_rotation, clump_ella, clump_ellb, clump_rotation, plot = 'off'):
    
    #Reproduce Maund/Reilly Model
    start = time.time() #record start time 
    
    #create grid of photons (identify each photon(pts) by x,y coordinate)
    nx, ny = nx, ny #1003, 1003 #number of photons(pts) in x and y directions, for a center pt nx = whole Number*xlim +1 (same for ny)
    xlim, ylim = 6, 6 #grid size (arbirary size (using units of AU) to fit photosphere)
    x = np.linspace(-xlim, xlim, nx) #values that grid will cover (using units of AU)
    y = np.linspace(-ylim, ylim, ny)
    centerx, centery = 0, 0 #specify center 
    
    #For photosphere ---------------------------------------------------------------------------------------------------------
    #plot ellipse photosphere radius
    ella, ellb = phot_ella, phot_ellb #5, 4.55 #ellipse radius ella = x, ellb = y so axis ratio: ellb/ella = 0.91 ratio from Reilly 2016 (arbitrary values - somewhat from homologous expansion caluclation)
    ell_thetas = np.linspace(0, 2*math.pi, 100) #angles - theta = 0 = +x , theta = pi = -x, theta = pi/2 = +y, theta = 3pi/2 = -y
    rotation = 2*math.pi-np.radians(phot_rotation) #photosphere major axis rotated 119degrees from Reilly 2016 rotates from x = 0 degrees (360-119 degrees but runs in radian)
    Ell = np.array([ella*np.cos(ell_thetas) , ellb*np.sin(ell_thetas)]) #ellipse radius pts in [x, y] array
    R_rot = np.array([[np.cos(rotation) , -np.sin(rotation)], [np.sin(rotation) , np.cos(rotation)]]) #2-D rotation matrix
    Ell_rot = np.zeros((2,Ell.shape[1])) #empty [x,y] array for rotated radius pts 
    for i in range(Ell.shape[1]):
        Ell_rot[:,i] = np.dot(R_rot,Ell[:,i]) #rotate ellipse by dot product of rotation matrix and original ellipse x,y coords array
    
    phot_xcoords = centerx+Ell_rot[0,:] #photosphere ellispe coordinates 
    phot_ycoords = centery+Ell_rot[1,:]
    
    #hold p, PA, I and x,y coords for each photon (x, y pt) in photosphere
    p_in_phot = []
    I_in_phot = []
    x_in_phot = []
    y_in_phot = []
    PA_in_phot = []
    zp = []; zI =[] #arrays for %p and I values for each pt in grid (z dimention for countour plot)
    #assign p values to pts in grid based on radius
    for rx in x:
        dummy = [] #dummy temp lists for countour plots 
        dum = []
        for ry in y: 
            r = np.sqrt(rx**2 + ry**2) #calculate radius at current pt
            theta = math.atan2(ry, rx) #calculate theta at current pt
            r_phot = np.sqrt((centerx+ella*np.cos(theta-rotation))**2 + (centery+ellb*np.sin(theta-rotation))**2) #calculate radius at pt on photosphere at same angle as pt in grid (account for rotated photosphere) 
            ###calc p here for contour plot otherwise saves time to do it only for pts in phot
            #p = 0.15*(r/r_phot); #calc p here for countour plot
            #I = 0.5*(1-np.sqrt(1-((r/r_phot)**2))); #calc I here for countour plot
            if r < r_phot:
                rphotx = r_phot*np.cos(theta); rphoty = r_phot*np.sin(theta) #calculate corresponding x,y coord on photosphere radius to rx,ry coord for pt (photon) in photosphere (r_phot already accounts for rotated photosphere)
                p = 0.15*(r/r_phot) #calculate %p equation from Maund 2010 (constant for %15 of light at limb polarized)
                I = 0.5*(1-np.sqrt(1-((r/r_phot)**2))) #I = Ir/I_0, the ratio of the lumonosity, using k = 0.5 for limb darkening from Maund 2010
                p_in_phot.append(p) #keep track of %p for pts in the photosphere
                I_in_phot.append(I) #keep track of Ir/I_0 for pts in the photosphere 
                x_in_phot.append(rx); y_in_phot.append(ry) #keep track of x, y coords. for pts in photosphere
                #calculate PA's used for model (16.5% parallel to photosphere, the rest randomly assigned)
                grn = random.uniform(0, 1) #generate random number to identify each pt (or photon) in phot by
                if grn > 0.165: # when random number > 0.165 assign random angle (between 0 and 180) to list
                    #grn2 = random.uniform(0, 1) #generate new random number to calculate random PA 
                    #PA_in_phot.append(180*grn2) # calculate and assign random angle (between 0 and 180) to list 
                    PA_in_phot.append(random.uniform(0, 180))#generate random angle between 0 and 180
                else:
                    #if random number < 0.165 append PA parallel to limb - this should cover that 16.5% of PA's be parallel to the limb Reilly 2016
                    #find slope for pts in ellipse according to https://www.anirdesh.com/math/algebra/ellipse-tangents.php
                    m = (ellb*rphotx)/(ella*np.sqrt(ella**2 - rphotx**2)) #calculate slope at corresponding location on ellipse 
                    if (rx <= 0 and ry >= 0) or (rx >= 0 and ry <= 0):
                        m = np.abs(m) #pts with positive slope (for ellipse rotated clockwise from major axis at x=0):
                    else: 
                        #negative slope
                        if m > 0: 
                            m = m*-1
                    angle = math.degrees(math.atan(m))#get angle of slope
                    if angle < 0:
                        angle = angle + 180 #add 180 to negative so in QU plane (0-180)
                    PA_in_phot.append(angle) #assigns every pt (photon) a PA that is the angle of slope in degrees
            #dummy.append(p) #(for countour plot)
            #dum.append(I)
        #zp.append(dummy) #(for countour plot)
        #zI.append(dum)
    
    #calculate q*I and u*I for pts in photosphere (do all math in intrinsic units)
    qI = (p_in_phot*np.cos(2*np.radians(PA_in_phot)))*I_in_phot  
    uI = (p_in_phot*np.sin(2*np.radians(PA_in_phot)))*I_in_phot
    
    #Create Table of values for pts (photons) in photosphere
    In_Phot = Table([x_in_phot, y_in_phot, p_in_phot, PA_in_phot, qI, uI, I_in_phot], names = ['x', 'y', 'p', 'pa', 'qI', 'uI', 'I'])
    
    #For obscuring clump ----------------------------------------------------------------------------------------------------
    ella_clump, ellb_clump = clump_ella, clump_ellb #3.5, 2.485 #ellipse radius ella = x, ellb = y so clump_ella = .88*phot_ella and ellb/ella = 0.71 from Reilly 2016
    ell_thetas_clump = np.linspace(0, 2*math.pi, 100) #angles - theta = 0 = +x , theta = pi = -x, theta = pi/2 = +y, theta = 3pi/2 = -y
    rotation_clump = 2*math.pi-np.radians(clump_rotation) #clump major axis rotated 160 degrees from Reilly 2016 rotates from x = 0 degrees 
    Ell_clump = np.array([ella_clump*np.cos(ell_thetas_clump) , ellb_clump*np.sin(ell_thetas_clump)]) #clump radius x,y coords array 
    R_rot_clump = np.array([[np.cos(rotation_clump) , -np.sin(rotation_clump)], [np.sin(rotation_clump) , np.cos(rotation_clump)]]) #2-D rotation matrix
    Ell_rot_clump = np.zeros((2, Ell_clump.shape[1])) #empty rotated clump radius x,y coords array 
    for i in range(Ell_clump.shape[1]):
        Ell_rot_clump[:,i] = np.dot(R_rot_clump, Ell_clump[:,i]) #rotate x,y coords and fill in empty array
    clump_xcoords = centerx+Ell_rot_clump[0,:] #rotated clump radius x coords, displaced from center 
    clump_ycoords = centery+Ell_rot_clump[1,:] #rotated clump radius x coords, displaced from center
    
    #Isolate data for pts in clump going through only photosphere pts (ie clump cannot extend past the bounds of the photosphere)
    In_Clump = Table(names = ['x', 'y', 'p', 'pa', 'qI', 'uI', 'I']) #creat table for pts within clump radius
    for ph_dat in In_Phot:
        r = np.sqrt(ph_dat['x']**2 + ph_dat['y']**2) #calculate radius at current pt in photosphere
        theta_c = math.atan2(ph_dat['y'], ph_dat['x']) - rotation_clump #calculate theta at current pt (account for clump rotation)
        r_clump = np.sqrt((centerx+ella_clump*np.cos(theta_c))**2 + (centery+ellb_clump*np.sin(theta_c))**2) #calculate radius of clump at same angle as current pt in photosphere
        if r < r_clump:
           In_Clump.add_row(ph_dat) #add values for pts inside clump radius into table
    
    #Isolate data in photosphere but not obscured by clump ------------------------------------------------------------------------
    Phot_noClump = setdiff(In_Phot, In_Clump)
    #calculate %p, PA in photosphere region not obscured by clump 
    #*NOTE* this value changes because the PA angle of each pt in grid is randomly assigned (except for the 16.5% that is parallel to the photosphere)
    Tot_P = (np.sqrt(np.sum(Phot_noClump['qI'])**2 + np.sum(Phot_noClump['uI'])**2))/np.sum(Phot_noClump['I']) #sum total q and u values by vector decomp. & divide by total I to get total p  
    qsum = np.sum(Phot_noClump['qI'])/np.sum(Phot_noClump['I']) #get q and u values to calculate new PAs
    usum = np.sum(Phot_noClump['uI'])/np.sum(Phot_noClump['I'])
    Tot_PA = np.degrees(.5*np.arctan2(usum, qsum)) #calculate PA by q and u vals 
    if Tot_PA < 0:
        Tot_PA = Tot_PA + 180
    #print("Total % p = " + str(Tot_P*100) + "Position Angle =" + str(Tot_PA)) #Values from Reilly are p = 0.7 +/- 0.5, theta = 70 +/- 3
    
    if plot == 'on':
        ####Plot model set up to check (only for few number pf photons (nx, ny)) 
        #make meshgrid from photon pts 
        xs, ys = np.meshgrid(x, y, sparse = True)#x, y values for countour plot
        xp, yp = np.meshgrid(x, y) #x, y values for individual grid pts
        fig, ax = plt.subplots(1)
        ax.set_xticks(np.arange(-xlim, xlim+1, 2)); ax.set_yticks(np.arange(-ylim, ylim+1, 2)); ax.axis('square')
        ax.plot(xp, yp, marker='o', color='k', linestyle='none') #plot grid pts (individual photons)
        ax.plot(centerx, centery, marker='o', color='r', linestyle='none') #plot center
        ax.plot(In_Phot['x'], In_Phot['y'], marker='o', color='m', linestyle='none') #plot pts (photons) in photosphere 
        ax.plot(In_Clump['x'], In_Clump['y'], marker='o', color='lime', linestyle='none') #plot pts (photons) in clump 
        ax.plot(Phot_noClump['x'], Phot_noClump['y'], marker='o', color='lightblue', linestyle='none') #plot pts (photons) in photosphere but not in clump 
        ax.plot(phot_xcoords, phot_ycoords, color = 'r') #plot photosphere radius
        ax.plot(clump_xcoords , clump_ycoords, 'y')#plot clump radius
        ###countour plot of I or P values by radius 
        #plot %p (or I) values as countour plot (one must be commented out at a time)
        #plot = ax.contourf(x, y, zp)#plot = ax.contourf(x, y, zI)
        #cbar = fig.colorbar(plot, orientation='horizontal'); cbar.set_label('% P')#cbar = fig.colorbar(plot, orientation='horizontal'); cbar.set_label('I')
    
    print('Model ran with ', nx*ny, 'total photons')
    print('Photosphere Ellipticity: ', phot_ellb/phot_ella, ', rotated by: ', phot_rotation, '°')
    print('Clump coverage: ', ((clump_ellb/clump_ella)/(phot_ellb/phot_ella))*100, '% of Photosphere')
    print('Clump Ellipticity: ', clump_ellb/clump_ella, ', rotated by: ', clump_rotation, '°')
    print('Returns: (Total P(%), Total PA(°), Total Photons)')
    end = time.time(); time_min = (end-start)/60; print("Model run time:", (time_min), "minutes")
    print('')
    
    return(Tot_P*100, Tot_PA, nx*ny)    
"""