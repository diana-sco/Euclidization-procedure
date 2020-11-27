#Tools to create, read, write image, do convolution, do deconvolution

import os
import galsim
import matplotlib.pyplot as plt
import math
import numpy as np

import config as config
import input as inp


def create_mock(x):
    """
    Creates a galaxy image for tests
    Returns a galsim image object
    """
    #image = galsim.ImageF(stampsize, stampsize)
    hlr = np.arange(inp.mock_reff_min, inp.mock_reff_max, inp.mock_reff_step)
    n = np.arange(inp.mock_n_min, inp.mock_n_max, inp.mock_n_step)
    a = x//len(hlr)
    b = x-a*len(hlr) 
    big_fft_params = galsim.GSParams(maximum_fft_size=90000)
    
    gal = galsim.Sersic(half_light_radius=hlr[b], n=n[a], flux=config.mock_flux, gsparams=big_fft_params)
    
    #gal = gal.shear(g1=g1, g2=g2)
    #gal.drawImage(image=image, scale=pixelscale, method='no_pixel')
    return gal, hlr[b], n[a]

def read_as_interpolatedimage(filepath, pixelscale, hdu=None):
    """
    Reads a FITS image into a galsim interpolatedImage object
    """
    image = galsim.fits.read(filepath, hdu=hdu)
    int_image = galsim.InterpolatedImage(image, scale=pixelscale)
    return int_image


def conv_with_PSF(image, psf):
    """
    image is a Galsim.image object
    psf is a Galsim.interpolatedImage object
    """
    convolut_image = galsim.Convolve([image, psf])
    return convolut_image


def deconv_with_PSF(image, psf):
    """
    image is a Galsim.image object
    psf is a Galsim.interpolatedImage object
    """
    inv_psf = galsim.Deconvolve(psf)
    deconv_image = galsim.Convolve(inv_psf, image)
    return deconv_image


def draw_image_conv_noise(convolut_image, aim, size_image):
    """
    Draw a FITS image with noise
    """
    if aim=='Eu_conv':
        gain=config.Eu_gain
        read_noise=config.Eu_rn
        expt=config.Eu_expt
        sky_level=config.Eu_sky
    elif aim=='HST_conv':
        gain=config.HST_gain
        read_noise=config.HST_rn
        expt=config.HST_expt
        sky_level=config.HST_sky
    shape = size_image.array.shape
    im = galsim.Image(shape[0], shape[1], scale=size_image.scale)
    im_empty = galsim.Image(shape[0], shape[1], scale=size_image.scale)
    convolut_image.drawImage(image=im)
    rng = galsim.BaseDeviate(config.random_seed)
    noise = galsim.CCDNoise(rng, sky_level=sky_level, gain=gain, read_noise=read_noise)
    im.addNoise(noise)
    im_empty.addNoise(noise)
    obj_noise = galsim.CorrelatedNoise(im_empty, rng=rng, scale=im.scale)
    return im, obj_noise


def symmetrize(image, obj_noise, size_image):
    shape = size_image.array.shape
    im = galsim.Image(shape[0], shape[1], scale=size_image.scale)
    image.drawImage(image=im)
    obj_noise.symmetrizeImage(im, order=4) #[default=4]
    return im
    
    
def draw_image_extranoise(image, size_image):
    """
    Draw a FITS image
    """
    shape = size_image.array.shape
    im = galsim.Image(shape[0], shape[1], scale=size_image.scale)
    im = image.drawImage(image=im)
    extra_noise_var = config.Eu_var - np.var(im.array[:5,:])
    rng = galsim.BaseDeviate(config.random_seed)
    gaussian_noise = galsim.GaussianNoise(rng, sigma=math.sqrt(extra_noise_var))
    im.addNoise(gaussian_noise)
    return im









    
    


    

