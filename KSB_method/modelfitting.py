''' Model fitting'''

from astropy.modeling import models, fitting
from astropy.stats import sigma_clip
from astropy.modeling import Fittable2DModel, Parameter
from astropy.io import fits
from scipy.special import gamma, gammaincinv
from scipy.stats import chisquare
import os
import math
import numpy as np
import config as config


class EllipSersic2D(Fittable2DModel):
	"""Custom Sersic Model, with ellipticity defined using (g1, g2), and total flux instead of amplitude
	
	This is based on astropy's Sersic2D from here:
	http://docs.astropy.org/en/stable/api/astropy.modeling.functional_models.Sersic2D.html#astropy.modeling.functional_models.Sersic2D
    
    Relations for the total flux:
    http://ned.ipac.caltech.edu/level5/March05/Graham/Graham2.html
    
    Instructions on how to build new models:
    http://docs.astropy.org/en/stable/modeling/new.html
	
	"""
	
	flux = Parameter(default=1.0, min=0.0) 
	r_eff = Parameter(default=1.0, min=0.0001, max=150.0)  #test also max=200 and max=150 in origin 100.0
	n = Parameter(default=4.0, min=0.1, max=6.0)
	x_0 = Parameter(default=0.0)
	y_0 = Parameter(default=0.0)
	g1 = Parameter(default=0.0, min=-1.0, max=1.0)
	g2 = Parameter(default=0.0, min=-1.0, max=1.0)
   
	@staticmethod
	def evaluate(x_array, y_array, flux, r_eff, n, x_0, y_0, g1, g2):
		
		theta = 0.5 * np.arctan2(g2, g1)
		g = np.hypot(g1, g2)
		#a, b = r_eff, r_eff*(1-g)/(1+g)  Original definition of Astropy from Malte
		a,b = r_eff/np.sqrt((1-g)/(1+g)), r_eff*np.sqrt((1-g)/(1+g)) #New definition for r_eff in order to have a match between mock_hlr and bestfit_reff
		
		bn = gammaincinv(2.0*n, 0.5)
		flux_n_factor =  n * np.exp(bn) * bn**(-2.0*n) * gamma(2.0*n)
		amplitude = flux / (2.0 * np.pi * flux_n_factor * a * b) # The amplitude from astropy's Sersic2D
		
		cos_theta, sin_theta = np.cos(theta), np.sin(theta)
		x_maj = (x_array - x_0) * cos_theta + (y_array - y_0) * sin_theta
		x_min = -(x_array - x_0) * sin_theta + (y_array - y_0) * cos_theta
		z = np.hypot(x_maj/a, x_min/b)
		
		return amplitude * np.exp(-bn * (z ** (1 / n) - 1))	
		
		
		
def make_a_fit(image, hlr, aim):
                if aim=='mock_image' or aim=='HST_decon':
                    gain=1.0
                    read_noise=1.0
                    expt=1.0
                elif aim=='Eu_conv' or aim=='Euclidiz':
                    gain=config.Eu_gain
                    read_noise=config.Eu_rn
                    expt=config.Eu_expt
                elif aim=='HST_conv':
                    gain=config.HST_gain
                    read_noise=config.HST_rn
                    expt=config.HST_expt


                #Make the fit
                stamp = image*1000 #here 10000 is meant to do the fit in any case, even if the flux is really low and so the failing of the fit
                weights = np.ones(stamp.shape)
                shape = image.shape
                x_array, y_array = np.meshgrid(np.arange(shape[0]), np.arange(shape[1]))
                init_flux = np.sum(stamp)
                
                ini_mod = EllipSersic2D(flux=init_flux, r_eff=config.init_reff, n=config.init_n, x_0=shape[0]/2, y_0=shape[1]/2, g1=config.init_g1, g2=config.init_g2)
                
                fitter = fitting.LevMarLSQFitter()
                
                fit_mod = fitter(ini_mod, x_array, y_array, stamp, weights=weights, maxiter = 1000, acc=1.0e-7, epsilon=1.0e-6, estimate_jacobian=False)
                
                
                #Do the calculation of the noise (sky level) for each image during the Euclidization procedure, considering a stripe of 1000 pixels with height 2 in the bottom part of the image. It is chosen enough far from the edge. This is used to calculate the chi2 and the reduced chi2
                var_image_chi = np.var(stamp[10:11,5:-5])
                
                
                ##Do the calculation of the noise (sky level) for each image during the Euclidization procedure, considering a stripe of 1000 pixels with height 2 in the bottom part of the image. It is chosen enough far from the edge. This is used to calculate the SNR in order to remove the *1000 we used in the stamp=image*1000
                var_image = np.var(image[10:11,5:-5])
                
        
                #Calculate the Chi2 over the pixels of the image, considering as variance the background sky level
                model = np.array(fit_mod(x_array, y_array).flatten())
                data = np.array(stamp.flatten())
                chi2= sum((data - model)**2/var_image_chi)
                red_chi2= chi2/config.truth_size**2
                
                #Calculate the SNR on the enter image of each step
                signal = (np.sum(image))*gain
                noises = (var_image*gain+(read_noise**2))
                aperture = math.pi*9*(fit_mod.r_eff[0]**2)
                
                if signal+aperture*noises > 0 and signal > 0:
                    SNR = signal/math.sqrt(signal+aperture*noises)
                else:
                    SNR = 0.0
                
                return fit_mod, fitter, chi2, red_chi2, SNR, init_flux
            
