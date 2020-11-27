#Script that creates dictionary to summarize the result from the fitting and KSB method estimations

from astropy.io import fits
from fitsio import FITS
import numpy as np
from astropy.table import Table
from astropy import units as u
import config as config
import math
import input as inp



def summary_result(aim, mock_n, reff, mock_g1, mock_g2, dx, dy, lens_g1, lens_g2, init_flux, KSB, ID):
    if aim=='mock_image' or aim=='HST_decon':
        gain=1.0
        zero_point=1.0
        expt=1.0
        read_noise=1.0
    elif aim=='Eu_conv' or aim=='Euclidiz':
        gain=config.Eu_gain
        zero_point=config.Eu_zp
        expt=config.Eu_expt
        read_noise=config.Eu_rn
    elif aim=='HST_conv':
        gain=config.HST_gain
        zero_point=config.HST_zp
        expt=config.HST_expt
        read_noise=config.HST_rn

    
    #Dictionary to store a row after each run of the code
    dict=np.zeros(1, dtype=[('aim', 'S32'), ('truth_pxscale', 'f8'), ('truth_size','f8'), ('mock_idx','f8'), ('mock_flux_mag','f8'), ('mock_SNR','f8'), ('mock_mu','f8'), ('mock_n','f8'), ('mock_reff','f8'), ('mock_g1','f8'), ('mock_g2','f8'), ('mock_flux','f8'), ('dx','f8'), ('dy','f8'), ('lens_g1','f8'), ('lens_g2','f8'), ('Eu_rn','f8'), ('Eu_gain','f8'), ('Eu_expt','f8'), ('HST_rn','f8'), ('HST_gain','f8'), ('HST_expt','f8'), ('init_n','f8'), ('init_reff','f8'), ('init_flux','f8'), ('init_flux_mag','f8'), ('g1_obs', 'f8'), ('g2_obs', 'f8'), ('g1_cor', 'f8'), ('g2_cor', 'f8'), ('shape_cor_err', 'f8'), ('m_sigma', 'f8'), ('m_amp', 'f8'), ('m_rho4', 'f8'), ('m_n', 'f8'), ('KSB_status', 'f8')])

    dict['aim']=aim
    dict['truth_pxscale']=config.truth_pxscale
    dict['truth_size']=config.truth_size
    dict['mock_idx']=inp.idx
    dict['mock_flux_mag']=inp.mock_mag
    dict['mock_SNR']=inp.Eu_snr
    dict['mock_mu']=config.mock_mu
    dict['mock_n']='%.2f'%mock_n
    dict['mock_reff']='%.2f'%reff
    dict['mock_g1']='%.4f'%mock_g1
    dict['mock_g2']='%.4f'%mock_g2
    dict['mock_flux']='%.1f'%config.mock_flux
    dict['dx']=dx
    dict['dy']=dy
    dict['lens_g1']='%.5f'%lens_g1
    dict['lens_g2']='%.5f'%lens_g2
    dict['Eu_rn']=config.Eu_rn
    dict['Eu_gain']=config.Eu_gain
    dict['Eu_expt']=config.Eu_expt
    dict['HST_rn']=config.HST_rn
    dict['HST_gain']=config.HST_gain
    dict['HST_expt']=config.HST_expt
    dict['init_n']=config.init_n
    dict['init_reff']=config.init_reff
    dict['init_flux']=init_flux/1000 # 10000 is the same value you can find in the model fitting to avoid the fitting failure
    dict['init_flux_mag']='%.2f'%(zero_point - 2.5*np.log10((init_flux/1000)*gain/expt))
    #dict['bestfit_flux']='%.5f'%(fit_mod.flux[0]/1000)  # 10000 is the same value you can find in the model fitting to avoid the fitting failure
    #dict['bestfit_flux_mag']='%.2f'%(zero_point - 2.5*np.log10((fit_mod.flux[0]/1000)*gain/expt))
    #dict['bestfit_reff_pix']='%.5f'%fit_mod.r_eff[0]
    #dict['bestfit_reff_arcsec']='%.5f'%(config.truth_pxscale*fit_mod.r_eff[0])   
    #dict['bestfit_n']='%.5f'%fit_mod.n[0]
    #dict['bestfit_x0']='%.5f'%fit_mod.x_0[0]
    #dict['bestfit_y0']='%.5f'%fit_mod.y_0[0]
    #dict['bestfit_g1']='%.5f'%fit_mod.g1[0]
    #dict['bestfit_g2']='%.5f'%fit_mod.g2[0]
    dict['g1_obs']='%.5f'%KSB.observed_shape.g1
    dict['g2_obs']='%.5f'%KSB.observed_shape.g2
    dict['g1_cor']='%.5f'%KSB.corrected_g1
    dict['g2_cor']='%.5f'%KSB.corrected_g2
    dict['shape_cor_err']='%.5f'%KSB.corrected_shape_err
    dict['m_sigma']='%.5f'%KSB.moments_sigma
    dict['m_amp']='%.5f'%KSB.moments_amp
    dict['m_rho4']='%.5f'%KSB.moments_rho4
    dict['m_n']=float(KSB.moments_n_iter)
    #dict['chi2']=chi2
    #dict['red_chi2']=red_chi2
    #dict['SNR']=SNR
    #dict['nfev']='%.5f'%fitter.fit_info['nfev']
    #dict['fvec']=fitter.fit_info['fvec']
    #dict['ierr']='%.5f'%fitter.fit_info['ierr']
    dict['KSB_status']=float(KSB.moments_status)
    
    #Creation of fits table containing one lines of results
    table = make_table_file()
    table.write(config.result_filename+'_'+str(aim)+'_'+'ID'+str(ID)+config.fits)
    f=FITS(config.result_filename+'_'+str(aim)+'_'+'ID'+str(ID)+config.fits, 'rw')   
    f[-1].append(dict)
    f.close()

    return 


def make_table_file():
	list=['aim', 'truth_pxscale', 'truth_size', 'mock_idx', 'mock_flux_mag', 'mock_SNR', 'mock_mu', 'mock_n', 'mock_reff', 'mock_g1', 'mock_g2', 'mock_flux', 'dx', 'dy', 'lens_g1', 'lens_g2', 'Eu_rn', 'Eu_gain', 'Eu_expt', 'HST_rn', 'HST_gain', 'HST_expt', 'init_n', 'init_reff', 'init_flux', 'init_flux_mag', 'g1_obs', 'g2_obs', 'g1_cor', 'g2_cor', 'shape_cor_err', 'm_sigma', 'm_amp', 'm_rho4', 'm_n', 'KSB_status']
	t = Table(names=list)
	t['aim']=t['aim'].astype('S32')
	return t

	
def print_result(Table):
	Table.write(config.result_filename)
	
	
	
