# -*- coding: utf-8 -*-
"""
Created on Tue Aug 16 15:43:00 2022

@author: sabri
"""

import numpy as np
import matplotlib.pyplot as plt

from get_data import *
from Binning_routines import *
from plot_codes import *
from astropy.table import Table
import random 

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

def mc_error_adjust(err_list, val_list, inst_err = 0.1):
    """ Uses a Monte Carlo method of generating a random number to determine 
    what value within a +/- error range to adjust an individual data point by.
    Takes in a list of error values and a corresponding list of point values ex. qerr, q. 
    Returns a list of error corrected point values. """
    from astropy.table import Table
    
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
    
    if ax == None: #set axis if not given 
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

#getting data
data2 = get_txtFITS("Epoch_2", "sn2012au_allcomb.qsum.txt", "sn2012au_allcomb.q.txt", "sn2012au_allcomb.qsig.txt", "sn2012au_allcomb.flx.txt", "sn2012au_allcomb.u.txt", "sn2012au_allcomb.usig.txt")
epoch2_20 = Bin_data(data2, 2.5, 20)
e2 = plot_line_velspace(epoch2_20, 5876, 15000, 10**-5.7, "Epoch2(Day26)")
print(e2)
#plotting
fig, ax = plt.subplots(1, figsize = (10, 10), sharex = True, sharey = True); fig.subplots_adjust(right = .85,  wspace=-.33, hspace=0.05)
fig.text(0.5, 0.08, '% Q Polarization', ha='center'); fig.text(0.05, 0.5, '% U Polarization', va='center', rotation='vertical') #set shared axis labels 
QU(e2) #function plots QU
get_loopy(e2['q'], e2['u'], e2['qerr'], e2['uerr'], ax) 
ax.legend()

print(e2['q'].type)

#Four farthest loop points:
qmin = get_min(e2, 'q'); qmin = e2[qmin[1]]
qmax = get_max(e2, 'q'); qmax = e2[qmax[1]]
umin = get_min(e2, 'u'); umin = e2[umin[1]]
umax = get_max(e2, 'u'); umax = e2[umax[1]]

sep0qu = [0.10, 0.68, 0.0125, 0.013]
sep75qu = [0.13, 0.16, 0.016, 0.019]
sep150qu = [0.57, -1.02, 0.025, 0.016]
sep225qu = [1.29, 1.72, 0.028, 0.025]
sep300qu = [1.40, 1.11, 0.015, 0.025]

qs = [0.10, 0.13, 0.57, 1.29, 1.40]
us= [0.68, 0.16, -1.02, 1.72, 1.11]
qerrs = [0.0125, 0.016, 0.025, 0.028, 0.015]
uerrs = [0.013, 0.019, 0.016, 0.025, 0.025]

ax.scatter(qs, us, marker = 'x', c ='k')
#get_loopy(qs, us, qerrs, uerrs, ax)

