#This code is the main procedure to Euclidiz HST images - Euclidization pipeline -

import galsim
import os
import math
import numpy as np
from sympy import *
import time

import config as config
from tools import create_mock, read_as_interpolatedimage, conv_with_PSF, deconv_with_PSF, draw_image_conv_noise, symmetrize, draw_image_extranoise
import input as inp
from modelfitting  import make_a_fit
from summarize_result import make_table_file, print_result, summary_result

import matplotlib.pyplot as plt
import random as rd

#Define the ellipticities for the mock galaxies which have a gaussian distribution
def mock_g():
    a=None
    b=None
    while (type(a)!=float or type(b)!=float):
        a=rd.gauss(config.mock_mean, config.mock_rms)
        b=rd.gauss(config.mock_mean, config.mock_rms)
        if type(a)==float and type(b)==float:
            if np.hypot(a, b)<1.0 and abs(a)<1.0 and abs(b)<1.0:
                return a, b
            else:
                a=None
                b=None
 
#Main function for the Euclidization procedure 
def main_main(x, gal, psf):
    print(x)
    a = x//(inp.NR*len(np.arange(inp.mock_reff_min, inp.mock_reff_max, inp.mock_reff_step))*len(np.arange(inp.mock_n_min, inp.mock_n_max, inp.mock_n_step))) #quoziente senza resta of the NR*number of created mock (lenght of hlr * lenght of n)
    mock_image = gal[0]
    mock_reff=float(gal[1])
    mock_n=float(gal[2])
    
    #Set shear parameters
    lens_g1=np.arange(inp.lens_g1_min, inp.lens_g1_max, inp.lens_g1_step)[a]
    lens_g2=0.0
    
    #Shape cancellation noise
    mock_g1_rd, mock_g2_rd = mock_g()   #generate the values
    mock_g1=[-1*mock_g1_rd, mock_g1_rd] # list of values for shape cancellation noise 
    mock_g2=[-1*mock_g2_rd, mock_g2_rd] #list of values for shape cancellation noise
    for ellip in range(2): #for the shape cancellation noise (gal and gal 90deg rotated)
        #ID to identify the fits files corresponding to each galaxy 
        ID = x+(ellip*(inp.NR*len(np.arange(inp.mock_reff_min, inp.mock_reff_max, inp.mock_reff_step))*len(np.arange(inp.mock_n_min, inp.mock_n_max, inp.mock_n_step))*len(np.arange(inp.lens_g1_min, inp.lens_g1_max, inp.lens_g1_step))))
        #Positional shift 
        dx, dy = np.random.uniform(-1*config.dx_shift, config.dy_shift), np.random.uniform(-1*config.dx_shift, config.dy_shift) 
        mock_im =  galsim.ImageF(config.truth_size, config.truth_size)
        mock_image = mock_image.shear(g1=mock_g1[ellip], g2=mock_g2[ellip])
        mock_image = mock_image.shift(dx, dy)
        mock_image.drawImage(image=mock_im, scale=config.truth_pxscale, method='no_pixel')
        
        #In order to fix the SNR or the magnitude
        if inp.idx==0:
            mag_image=inp.mock_mag
        elif inp.idx==1:
            y= symbols('y')
            Eqflux= solveset(Eq(y*config.Eu_gain/(sqrt(y*config.Eu_gain+math.pi*9*((mock_reff/config.truth_pxscale)**2)*(config.Eu_var*config.Eu_gain+config.Eu_rn**2))), inp.Eu_snr), y)
            flux_Eu = float(Eqflux.args[0])
            mag_image = (config.Eu_zp - 2.5*np.log10((flux_Eu)*config.Eu_gain/config.Eu_expt))
        else:
            print('CHOOSE THE RIGHT INDEX!')

        flux_Eu=(config.Eu_expt/config.Eu_gain)*10**(-0.4*(mag_image-config.Eu_zp)) #ADU
        flux_HST=(config.HST_expt/config.HST_gain)*10**(-0.4*(mag_image-config.HST_zp)) #ADU
        
        #Interpolate the mock image
        mock_image=galsim.InterpolatedImage(mock_im, scale=config.truth_pxscale)
        
        #Convolution mock image with Euclid PSF and model fitting of convoluted image. I also add cosmic shear before convolution and CCD noise after the convolution. KSB method for shear measurements is applied
        mock_lens = mock_image.lens(lens_g1, lens_g2, config.lens_mu)  #add this where mu is 1 or 1.5 
        Eu_conv_image = conv_with_PSF(mock_lens, psf['Euclid_psf'])
        Eu_conv_image = Eu_conv_image*(flux_Eu/(Eu_conv_image.flux))
        Eu_conv_im, Eu_conv_noise = draw_image_conv_noise(Eu_conv_image, aim='Eu_conv', size_image=mock_im)
        #fit_mod_Eu, fitter_Eu, chi2_Eu, red_chi2_Eu, SNR_Eu, init_flux_Eu = make_a_fit(Eu_conv_im.array, mock_reff, aim='Eu_conv')
        init_flux_Eu = np.sum(Eu_conv_im.array)
        KSB_Eu = galsim.hsm.EstimateShear(Eu_conv_im, psf['im_Euclid_psf'], sky_var=config.Eu_var, shear_est='KSB', strict='False')
        summary_result(aim='Eu_conv',
                               mock_n=mock_n,
                               reff=mock_reff,
                               mock_g1=mock_g1[ellip],
                               mock_g2=mock_g2[ellip],
                               dx=dx,
                               dy=dy,
                               lens_g1=lens_g1,
                               lens_g2=lens_g2,
                               init_flux=init_flux_Eu,
                               #fit_mod=fit_mod_Eu,
                               #fitter=fitter_Eu,
                               KSB=KSB_Eu,
                               #chi2=chi2_Eu,
                               #red_chi2=red_chi2_Eu,
                               #SNR=SNR_Eu,
                               ID=ID)
        
        #Convolution mock image with HST PSF and model fitting of convoluted image. I also add CCD noise after the convolution. KSB method for shear measurements is applied
        HST_conv_image = conv_with_PSF(mock_image, psf['Hubble_psf'])
        HST_conv_image = HST_conv_image*(flux_HST/(HST_conv_image.flux))
        HST_conv_im, HST_conv_noise = draw_image_conv_noise(HST_conv_image, aim='HST_conv', size_image=mock_im)
        #fit_mod_HST, fitter_HST, chi2_HST, red_chi2_HST, SNR_HST, init_flux_HST = make_a_fit(HST_conv_im.array, mock_reff, aim='HST_conv')
        #KSB_HST = galsim.hsm.EstimateShear(HST_conv_im, psf['im_Hubble_psf'], sky_var=config.HST_var, shear_est='KSB', strict='False')
        #summary_result(aim='HST_conv',
                               #mock_n=mock_n,
                               #reff=mock_reff,
                               #mock_g1=mock_g1[ellip],                                               #mock_g2=mock_g2[ellip],
                               #dx=dx,
                               #dy=dy,
                               #lens_g1=lens_g1,
                               #lens_g2=lens_g2,
                               #init_flux=init_flux_HST,
                               #fit_mod=fit_mod_HST,
                               #fitter=fitter_HST,
                               #KSB=KSB_HST,
                               #chi2=chi2_HST,
                               #red_chi2=red_chi2_HST,
                               #SNR=SNR_HST,
                               #ID=ID)
        
        #Deconvolution of convoluted image with HST PSF 
        HST_image = galsim.InterpolatedImage(HST_conv_im, scale=HST_conv_im.scale)
        HST_deconv_image = deconv_with_PSF(HST_image, psf['Hubble_psf'])
        
        #Convolution of deconvoluted (HST IMAGE * HST PSF) image with Euclid PSF
        HST_deconv_lens = HST_deconv_image.lens(lens_g1, lens_g2, config.lens_mu) #add this mu=1.5 or 1 
        Euclidiz_image = conv_with_PSF(HST_deconv_lens, psf['Euclid_psf'])  
        
        #Symmetrize the Euclidiz image after the deconvolution by Euclid PSF and rescale for the flux. I add also some extra gaussian noise. Fittong and KSB method for shear measurements are applied
        Euclidiz_sym_im = symmetrize(Euclidiz_image, HST_conv_noise, size_image=HST_conv_im)
        Euclidiz = galsim.InterpolatedImage(Euclidiz_sym_im, scale=Euclidiz_sym_im.scale)
        Euclidiz = Euclidiz*(flux_Eu/flux_HST)
        Euclidiz_im = draw_image_extranoise(Euclidiz, size_image=Euclidiz_sym_im)
        #fit_mod_Euz, fitter_Euz, chi2_Euz, red_chi2_Euz, SNR_Euz, init_flux_Euz = make_a_fit(Euclidiz_im.array, mock_reff, aim='Euclidiz')
        init_flux_Euz = np.sum(Euclidiz_im.array) 
        KSB_Euz = galsim.hsm.EstimateShear(Euclidiz_im, psf['im_Euclid_psf'], sky_var=config.Eu_var, shear_est='KSB', strict='False')
        summary_result(aim='Euclidiz',
                       mock_n=mock_n,
                               reff=mock_reff,
                               mock_g1=mock_g1[ellip],
                               mock_g2=mock_g2[ellip],
                               dx=dx,
                               dy=dy,
                               lens_g1=lens_g1,
                               lens_g2=lens_g2,
                               init_flux=init_flux_Euz,
                               #fit_mod=fit_mod_Euz,
                               #fitter=fitter_Euz,
                               KSB=KSB_Euz,
                               #chi2=chi2_Euz,
                               #red_chi2=red_chi2_Euz,
                               #SNR=SNR_Euz,
                               ID=ID)
    return
                           



    
    
    
