#Input parameters for the Euclidization procedure that you can change for each run

#mock parameters
mock_reff_min=0.2   #arcsec  
mock_reff_max= 0.41  #arcsec  
mock_reff_step=0.2  #arcsec

mock_n_min=1.0
mock_n_max= 3.01 #4.5
mock_n_step=2.0

#shear parameters
lens_g1_min=-0.04 #-0.1
lens_g1_max=0.041
lens_g1_step= 0.004 #0.004

lens_g2_min=0.1 #-0.1
lens_g2_max=0.1001
lens_g2_step=0.004

#Number realization=NR
NR=14000

#For images with fixed magnitude of the mock galaxy, set idx=0; for images with fixed SNR, set idx=1
idx= 1 
mock_mag=20.0 #mag
Eu_snr=20

#This is the directory containing the fits files in output
workdir= 'SNR20_3t_NR14000_g1' 
