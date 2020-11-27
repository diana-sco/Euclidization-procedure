#Code to stack the output fit tables (result) in one fits table
import os
from astropy.io import fits
from astropy.table import Table, vstack
import multiprocessing as mlp
from time import perf_counter

import input as inp

def read_fits(x, filename):
    table = Table.read(inp.workdir+'/'+str(filename), format='fits')
    return table
    
    
#Time counter
t0=perf_counter()
a = os.listdir(inp.workdir)
print('1')
#Pool to read all the fits tables
pool = mlp.Pool(mlp.cpu_count())
t = pool.starmap(read_fits, [(x, a[x]) for x in range(len(a))])
pool.close()
pool.join()
print('2')
table = vstack(t)
print('3')
table.write('prova.fits')

t1=perf_counter() 
print('************'+str(t1-t0))

