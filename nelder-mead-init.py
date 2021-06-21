#!/usr/bin/env python3

import numpy as np
import copy
import matscipy.calculators.eam.calculator
import shelve
import os

def find_peak():
    os.system('gen-ns')
    os.system('./ns_analyse 64Cu.energies -M 0 -n 600 -D 5 > pf')
    T,Cp = numpy.loadtxt('pf', usecols=(0,3), unpack=True)
    index_PT =  numpy.where(Cp == numpy.amax(Cp))
    print(int(T[index_PT]))
    return(experimental - int(T[index_PT]))

def gen_poten(param):
    potential = matscipy.calculators.eam.io.read_eam('Cu_Zhou04.eam.alloy', 'eam/alloy')
    matscipy.calculators.eam.io.write_eam(potential[0], potential[1], (potential[2] * param[0]), (potential[3] * param[1]), (potential[4] *  param[2]), 'test.eam.alloy','eam/alloy')

def save_variables():
    cwd = os.getcwd()
    filename= str(cwd) + '/variables'
    shelf = shelve.open(filename, 'n')
    for key in globals():
        try:
            shelf[key] = globals()[key]
        except:
            print('ERROR shleving: {0}'.format(key))
    shelf.close()

def load_variables():
    cwd = os.(getcwd)
    shelf = shelf.open(str(cwd +'/variables'))
    for key in shelf:
        globals()[key] = shelf[key]
    shelf.close

# Read in and modify initial potential, write new potential files




#Perform nested sampling optimisation, we want to minimise score which is the difference of PT from experimental value.


step = 0.05
no_improv_thr = 10e-6
n_improv_break= 10
max_iter = 1
alpha = 1.
gamma = 2.
rho = -0.5
sigma = 0.5

#my vars

x_start = np.array([1.0, 1.0, 1.0])
experimental = 1400

# init
dim = len(x_start)
save_variables()
gen_poten(x_start)
load_variables()
prev_best = find_peak(x_start)
no_improv = 0
res = [[x_start, prev_best]]

#for i in range(dim):
i=0
while i < len(dim)
    x = copy.copy(x_start)
    x[i] = x[i] + step
    save_variables()
    gen_poten(x)
    load_variables()
    score = find_peak(x)
    res.append([x, score])
    i += 1

# simplex iter
iters = 0
while 1:
    # order
    res.sort(key=lambda x: x[1])
    best = res[0][1]

    # break after max_iter
    if max_iter and iters >= max_iter:
        return res[0]
    iters += 1

    # break after no_improv_break iterations with no improvement
    print '...best so far:', best

    if best < prev_best - no_improve_thr:
        no_improv = 0
        prev_best = best
    else:
        no_improv += 1

    if no_improv >= no_improv_break:
        return res[0]

    # centroid
    x0 = [0.] * dim
    for tup in res[:-1]:
        for i, c in enumerate(tup[0]):
            x0[i] += c / (len(res)-1)

    # reflection
    xr = x0 + alpha*(x0 - res[-1][0])
    save_variables()
    gen_poten(xr)
    load_variables()
    rscore = find_peak(xr)
    if res[0][1] <= rscore < res[-2][1]:
        del res[-1]
        res.append([xr, rscore])
        continue

    # expansion
    if rscore < res[0][1]:
        xe = x0 + gamma*(x0 - res[-1][0])
        save_variables()
        gen_poten(xe)
        load_variables()
        escore = find_peak(xe)
        if escore < rscore:
            del res[-1]
            res.append([xe, escore])
            continue
        else:
            del res[-1]
            res.append([xr, rscore])
            continue

    # contraction
    xc = x0 + rho*(x0 - res[-1][0])
    save_variables()
    gen_poten(xc)
    load_variables()
    cscore = find_peak(xc)
    if cscore < res[-1][1]:
        del res[-1]
        res.append([xc, cscore])
        continue

    # reduction
    x1 = res[0][0]
    nres = []

    print('got to reduction, cant restart in a for loop')
    print('res:')
    print(res)
    exit()
    #for tup in res:
    i=0
    while i < len(res):
        redx = x1 + sigma*(tup[0] - x1)
        save_variables()
        gen_poten(redx)
        load_variables()
        score = find_peak(redx)
        nres.append([redx, score])
        i += 1
    res = nres






#def reload_variables():
#    variables = np.load('variables.npz')
#    res =       variables['res']
#    best =      variables['best']
#    iters =     variables['iters']
#    no_improv = variables['no_improv']
#    dim =       variables['dim']
#    x0 =        variables['x0'] 
#    xr =        variables['xr']
#    rscore =    variables['rscore']
#def save_variables():
#    np.savez('variables', res=res, best=best, iters=iters, no_improv=no_improv, dim=dim, x0=x0, xr=xr, rscore=rscore, xe=xe)
