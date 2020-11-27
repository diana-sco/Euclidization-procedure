#Configuration file containing all info for the Euclidization procedure scripts

import os
import input as inp
import numpy as np

#extension fits filename
fits='.fits'
png='.png'

#this is the directory containing the fits files in output
workdir=inp.workdir
if not os.path.isdir(workdir):
    os.mkdir(workdir)  
    
#this is the directory containing the png files in output
workdir_png= workdir+'_png'
if not os.path.isdir(workdir_png):
    os.mkdir(workdir_png)   

#values for the drawn images (fits file in output)
truth_pxscale=0.03  #arcsec/pixels
truth_size=512 #pixels 512 x 512

#Seed to create noise
#random_seed = 1314662
random_seed = None #in order to have different realization of noise

#Noise and SNR --- Euclid set-up
Eu_pxs = 0.1 #arcsec/pixels
Eu_gain = 3.1 #1.7 or 1000.0 electrons/ADU
Eu_rn = (4.2*truth_pxscale/Eu_pxs)  # readnoise rescale for the pixscale I am using (0.03), 4.2 isvalue from Malte's paper, old value= 0.3 electrons/pixels
Eu_expt=1695.0 #1695 or 565.0 seconds **TO CHANGE**
Eu_sky_mag=22.35 #mag/arcsec^2
Eu_zp=24.6 #mag
Eu_sky=(((Eu_expt/Eu_gain)*10**(-0.4*(Eu_sky_mag-Eu_zp)))*(truth_pxscale**2))  #ADU/pixel
Eu_var=(Eu_sky/Eu_gain)+(Eu_rn/Eu_gain)**2 #from Galsim documentation

#Noise and SNR --- HST set-up
HST_pxs = 0.05 #arcsec/pixels
HST_gain = 2.0 #1.7 or 1000.0 electrons/ADU
HST_rn = (5.0*truth_pxscale/HST_pxs)  #readnoise rescaled for the pixscale I am using (0.03),5.0 isvalue from documentation (it should be 2/3 times bigger than the Euclid's one)
HST_expt=1000.0 #seconds
HST_sky_mag=22.5 #mag that is the avarage btw 22 and 23 (values Tim said), 23.34 comes from the formula
HST_zp=26.5 #mag
HST_sky=(((HST_expt/HST_gain)*10**(-0.4*(HST_sky_mag-HST_zp)))*(truth_pxscale**2))  #ADU/pixel
HST_var=(HST_sky/HST_gain)+(HST_rn/HST_gain)**2 #from Galsim documentation

#mock parameters
mock_flux=10000  #this value is not involved righ now, it could be also 1.e4  100 or 2.e3 counts
mock_mu=1.0
#to create the gaussian distribution for the mock_g1 and mock_g2:
mock_mean=0.0
mock_rms=0.3

#'lens parameters
lens_mu=1.0

#Positional shift
dx_shift=0.015 #arcsec
dy_shift=0.015 #arcsec
    
#Euclid PSF
Eu_psf='sed_true_26892756.os.fits'
Eu_psf_pxs = 0.02 #arcsec/pixels
Eu_psf_size = 800 #pixels

'''HST PSF'''
HST_psf='tt3_64_64_1.5.fits'
HST_psf_pxs = 0.01    #arcsec/pixels
HST_psf_size = 251 #pixels
#fitter parameters
init_n=3.0
init_reff=5.0
init_x0=truth_size/2.0
init_y0=truth_size/2.0
init_g1=0.2
init_g2=0.1

#fits file output (table) cointaing the results
result_filename = os.path.join(workdir, workdir)
