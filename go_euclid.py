#this code parallelizes the procedure

import galsim
import multiprocessing as mlp
from time import perf_counter
from astropy.table import Table, vstack
import numpy as np
from functools import partial
from time import perf_counter
from fitsio import FITS

import input as inp
from summarize_result import make_table_file, print_result
from main import main_main
from tools import create_mock, read_as_interpolatedimage
import config as config



def main_func():
    #Time counter
    t0=perf_counter() 
 
    #Definition of the pool to create the mock galaxy as GSObject
    pool = mlp.Pool(mlp.cpu_count())
    gal = pool.map(create_mock, range(len(np.arange(inp.mock_reff_min, inp.mock_reff_max, inp.mock_reff_step))*
                                      len(np.arange(inp.mock_n_min, inp.mock_n_max, inp.mock_n_step))))
    pool.close()
    pool.join()
    print('************end creation mock***********')
    
    #Import and interpolate the PSFs
    Euclid_psf=read_as_interpolatedimage(config.Eu_psf, config.Eu_psf_pxs, hdu=1)
    im_Euclid_psf = galsim.Image(config.Eu_psf_size, config.Eu_psf_size, scale=config.Eu_psf_pxs)
    Euclid_psf.drawImage(image=im_Euclid_psf)
    Hubble_psf=read_as_interpolatedimage(config.HST_psf, config.HST_psf_pxs)
    im_Hubble_psf = galsim.Image(config.HST_psf_size, config.HST_psf_size, scale=config.HST_psf_pxs)
    Hubble_psf.drawImage(image=im_Hubble_psf)
    
    #Creation of a dictionary to have the PSFs as keywords in the main function
    psf = {}
    psf['Euclid_psf'] = Euclid_psf
    psf['im_Euclid_psf'] = im_Euclid_psf
    psf['Hubble_psf'] = Hubble_psf
    psf['im_Hubble_psf'] = im_Hubble_psf
    
    #Pool for the main_main function
    pool = mlp.Pool(mlp.cpu_count())
    pool.starmap(main_main, [(x, gal[x//(inp.NR*len(np.arange(inp.lens_g1_min, inp.lens_g1_max, inp.lens_g1_step)))], psf) for x in range(len(gal)*inp.NR*len(np.arange(inp.lens_g1_min, inp.lens_g1_max, inp.lens_g1_step)))]) 
    pool.close()
    pool.join()
    
    t1=perf_counter() 
    print('************'+str(t1-t0))
    return

if __name__=='__main__':
    main_func()
    

