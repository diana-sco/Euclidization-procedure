#Code to stack the output fit tables (result) in one fits table and do plots
import matplotlib.pyplot as plt
import numpy as np
import csv
import os
from astropy.io import fits
from astropy.table import Table, vstack
import multiprocessing as mlp
from time import perf_counter
import math
import statsmodels.api as sm

import input as inp

#Definiton to read ll fits tables
def read_fits(x, filename):
    #table = fits.open(inp.workdir+'/'+str(filename), format='fits')[1].data
    table = fits.open('prova/'+str(filename), format='fits')[1].data
    print(x)
    return table

    
#Time counter
t0=perf_counter()

#a = os.listdir(inp.workdir)
a = os.listdir('prova') #(UNCOMMENT IF YOU AHVE MORE THAN ONE FILES)
print('start')

#Pool to read all the fits tables (UNCOMMENT IF YOU AHVE MORE THAN ONE FILES)
pool = mlp.Pool(mlp.cpu_count())
data = pool.starmap(read_fits, [(x, a[x]) for x in range(len(a))]) 
pool.close()
pool.join()
print('end reading files')

#If you have only a fits table with all result inside it
#data = fits.open('prova_g1.fits')[1].data
#print('end reading files')

aim=['Eu_conv', 'Euclidiz']
lens_g1=sorted(list(set([float(data[i]['lens_g1']) for i in range(len(data))])))
lens_g2=sorted(list(set([float(data[i]['lens_g2']) for i in range(len(data))])))
mock_reff=sorted(list(set([float(data[i]['mock_reff']) for i in range(len(data))])))
mock_n=sorted(list(set([float(data[i]['mock_n']) for i in range(len(data))])))

N=100 #number of valid realizations (so, galaxies) for each datapoint

#Define global lists
euconvKSBmean=[]
euconvKSBstd=[]
euclidizKSBmean=[]
euclidizKSBstd=[]

#Pool to create a global list
for par_aim in aim:
    for par_lens_g1 in lens_g1:
        condition = [par_aim, par_lens_g1]
        lista = [float(d['g1_cor']) for d in data
                 if d['aim']==str(condition[0]) and float(d['lens_g1'])==float(condition[1]) and float(d['KSB_status'])==0.0 and np.hypot(float(d['g1_cor']), float(d['g2_cor']))<1]
        if par_aim=='Eu_conv':
            euconvKSBmean.append(np.mean(list(lista[0:(N*len(mock_n)*len(mock_reff))])))
            euconvKSBstd.append(np.std(list(lista[0:(N*len(mock_n)*len(mock_reff))]))/math.sqrt(N*len(mock_n)*len(mock_reff)))
        elif par_aim=='Euclidiz':
            euclidizKSBmean.append(np.mean(list(lista[0:(N*len(mock_n)*len(mock_reff))])))
            euclidizKSBstd.append(np.std(list(lista[0:(N*len(mock_n)*len(mock_reff))]))/math.sqrt(N*len(mock_n)*len(mock_reff))) 
print('end global list')  

#Define lists considering galaxy param          
par_euconvKSBmean=[]
par_euconvKSBstd=[]
par_euclidizKSBmean=[]
par_euclidizKSBstd=[]

for par_mock_n in mock_n:
    for par_mock_reff in mock_reff:
        par_euconvKSB_mean=[]
        par_euconvKSB_std=[]
        par_euclidizKSB_mean=[]
        par_euclidizKSB_std=[]
        for par_aim in aim:
            for par_lens_g1 in lens_g1:
                par_condition = [par_aim, par_lens_g1, par_mock_reff, par_mock_n]
                par_lista = [float(d['g1_cor']) for d in data
                             if d['aim']==str(par_condition[0]) and float(d['lens_g1'])==float(par_condition[1]) and float(d['mock_reff'])==float(par_condition[2]) and float(d['mock_n'])==float(par_condition[3]) and float(d['KSB_status'])==0.0 and np.hypot(float(d['g1_cor']),float(d['g2_cor']))<1]
                mean=np.mean(list(par_lista[0:N])) 
                std=np.std(list(par_lista[0:N]))/math.sqrt(N)
                print(len(par_lista))
                if par_aim=='Eu_conv':
                    par_euconvKSB_mean.append(mean)
                    par_euconvKSB_std.append(std)
                elif par_aim=='Euclidiz':
                    par_euclidizKSB_mean.append(mean)
                    par_euclidizKSB_std.append(std)
        par_euconvKSBmean.append([par_mock_n, par_mock_reff, par_euclidizKSB_mean])
        par_euconvKSBstd.append([par_mock_n, par_mock_reff, par_euclidizKSB_mean])
        par_euclidizKSBmean.append([par_mock_n, par_mock_reff, par_euclidizKSB_mean])
        par_euclidizKSBstd.append([par_mock_n, par_mock_reff, par_euclidizKSB_mean])
print('end par_list')
#***********************To plot quantities WHITHOUT any distinction of mock_reff and mock_n****************************
#To have a shift in the data point for visualization purpose        
lens_g1_plus=[]
lens_g1_minus=[]
for i in range(len(lens_g1)):
    lens_g1_plus.append(lens_g1[i]+0.0005)
    lens_g1_minus.append(lens_g1[i]-0.0005)
    
#Difference between KSB estimated and input g for Eu_conv and Euclidiz   
diff_Eu_conv=[]
for i in range(len(lens_g1)):
    diff_Eu_conv.append(euconvKSBmean[i] - lens_g1[i])
    
diff_Euz=[]
for i in range(len(lens_g1)):
    diff_Euz.append(euclidizKSBmean[i] - lens_g1[i])    
    
# To estimate m and c (biases) with errors
fit_err = sm.add_constant(lens_g1) #in the bracket the x variable
model_Euc = sm.OLS(diff_Eu_conv, fit_err) #in the bracket y, x
model_Euz = sm.OLS(diff_Euz, fit_err) #in the bracket y, x
results_Euc = model_Euc.fit()
results_Euz = model_Euz.fit()
print(results_Euc) 

#Linear fits for Eu_conv and Euclidiz
fit_line_Eu_conv=[]
for i in range(len(lens_g1)):
    fit_line_Eu_conv.append(results_Euc.params[1]*lens_g1[i]+results_Euc.params[0])
fit_line_Euz=[]
for i in range(len(lens_g1)):
    fit_line_Euz.append(results_Euz.params[1]*lens_g1[i]+results_Euz.params[0])    
        
#To do the plots
plt.errorbar(lens_g1_minus, diff_Eu_conv, euconvKSBstd, linestyle='', marker='o', markersize=5, color= 'r', label='Eu_conv'+'  m='+str('%.3f'%results_Euc.params[1])+'$\pm$'+str('%.3f'%results_Euc.bse[1])+'  c='+str('%.3f'%results_Euc.params[0])+'$\pm$'+str('%.3f'%results_Euc.bse[1]))
plt.plot(lens_g1, fit_line_Eu_conv, '-r')
plt.errorbar(lens_g1_plus, diff_Euz, euclidizKSBstd, linestyle='', marker='o', markersize=5, color= 'g', label='Euclidiz'+'  m='+str('%.3f'%results_Euz.params[1])+'$\pm$'+str('%.3f'%results_Euz.bse[1])+'  c='+str('%.3f'%results_Euz.params[0])+'$\pm$'+str('%.3f'%results_Euz.bse[1]))
plt.plot(lens_g1, fit_line_Euz, '-g')
plt.xlabel('g1$_{inp}$')   #TO CHANGE lens_g1
plt.ylabel('g1$_{KSB}$ - g1$_{inp}$')  #TO CHANGE lens_g1
plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
#plt.xlim(0.01, 0.05)
plt.legend()
plt.savefig('g1KSB_NR10000_3t_SNR20_glob.png')   #TO CHANGE lens_g1
plt.close()            
 
#***********************To plot quantities WITH DISTONCTION of mock_reff and mock_n**************************** 
bias = []
for j in range(len(par_euconvKSBmean)):
    d = {}
    #Difference between observed and input g for Eu_conv and Euclidiz   
    diff_Eu_conv=[]
    for i in range(len(lens_g1)):
        diff_Eu_conv.append(par_euconvKSBmean[j][2][i] - lens_g1[i])
    
    diff_Euz=[]
    for i in range(len(lens_g1)):
        diff_Euz.append(par_euclidizKSBmean[j][2][i] - lens_g1[i])      
    
    # To estimate m and c with errors
    fit_err = sm.add_constant(lens_g1) #in the bracket the x variable
    model_Euc = sm.OLS(diff_Eu_conv, fit_err) #in the bracket y, x
    model_Euz = sm.OLS(diff_Euz, fit_err) #in the bracket y, x
    results_Euc = model_Euc.fit()
    results_Euz = model_Euz.fit()
    
    #Linear fits for Eu_conv and Euclidiz
    fit_line_Eu_conv=[]
    for i in range(len(lens_g1)):
        fit_line_Eu_conv.append(results_Euc.params[1]*lens_g1[i]+results_Euc.params[0])
    fit_line_Euz=[]
    for i in range(len(lens_g1)):
        fit_line_Euz.append(results_Euz.params[1]*lens_g1[i]+results_Euz.params[0])  
        
    #Dictionary to store the multiplicative and additive bias
    d['m_Euc']=results_Euc.params[1]
    d['m_err_Euc']=results_Euc.bse[1]
    d['c_Euc']=results_Euc.params[0]
    d['c_err_Euc']=results_Euc.bse[1]
    d['m_Euz']=results_Euz.params[1]
    d['m_err_Euz']=results_Euz.bse[1]
    d['c_Euz']=results_Euz.params[0]
    d['c_err_Euz']=results_Euz.bse[1]
    d['mock_reff']=par_euconvKSBmean[j][1]
    d['mock_n']=par_euconvKSBmean[j][0]
    bias.append(d)

    #To do the plots
    plt.errorbar(lens_g1_minus, diff_Eu_conv, par_euconvKSBstd[j][2], linestyle='', marker='o', markersize=5, color= 'r', label='Eu_conv'+'  m='+str('%.3f'%results_Euc.params[1])+'$\pm$'+str('%.3f'%results_Euc.bse[1])+'  c='+str('%.3f'%results_Euc.params[0])+'$\pm$'+str('%.3f'%results_Euc.bse[1]))
    plt.plot(lens_g1, fit_line_Eu_conv, '-r')
    plt.errorbar(lens_g1_plus, diff_Euz, par_euclidizKSBstd[j][2], linestyle='', marker='o', markersize=5, color= 'g', label='Euclidiz'+'  m='+str('%.3f'%results_Euz.params[1])+'$\pm$'+str('%.3f'%results_Euz.bse[1])+'  c='+str('%.3f'%results_Euz.params[0])+'$\pm$'+str('%.3f'%results_Euz.bse[1]))
    plt.plot(lens_g1, fit_line_Euz, '-g')
    plt.title('mock_n: '+str(par_euconvKSBmean[j][0])+', mock_reff: '+str(par_euconvKSBmean[j][1]))
    plt.xlabel('g1$_{inp}$')   #TO CHANGE lens_g1
    plt.ylabel('g1$_{KSB}$ - g1$_{inp}$')  #TO CHANGE lens_g1
    plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
    #plt.xlim(0.01, 0.05)
    plt.legend()
    plt.savefig('g1KSB_NR10000_3t_SNR20_mockn'+str(par_euconvKSBmean[j][0])+'_mockreff'+str(par_euconvKSBmean[j][1])+'.png')   #TO CHANGE lens_g1
    plt.close()            
 
#******************To plot m_Euz-m_Euc (c_Euz-c_Euc) vs mock_reff and  mock_n******************************
#plot difference between multiplicative bias
for n in mock_n:
    plt.errorbar(mock_reff,
                 [(b['m_Euc']-b['m_Euz']) for b in bias if b['mock_n']==n],
                 [(b['m_err_Euc']+b['m_err_Euz']) for b in bias if b['mock_n']==n],
                 linestyle='--', marker='o', markersize=5, color= 'k',
                 label='mock_n = '+str(n))
plt.xlabel('mock_reff')
plt.ylabel('m_Euc-m_Euz')
plt.legend()
plt.savefig('m_diff_g1.png')        
plt.close()            
    
    
#plot difference between additive bias
for n in mock_n:
    plt.errorbar(mock_reff,
                 [(b['c_Euc']-b['c_Euz']) for b in bias if b['mock_n']==n],
                 [(b['c_err_Euc']+b['c_err_Euz']) for b in bias if b['mock_n']==n],
                 linestyle='--', marker='o', markersize=5, color= 'b',
                 label='mock_n = '+str(n))
plt.xlabel('mock_reff')
plt.ylabel('c_Euc-c_Euz')
plt.legend()
plt.savefig('c_diff_g1.png')        
plt.close()            
    
t1=perf_counter() 
print('************'+str(t1-t0))

