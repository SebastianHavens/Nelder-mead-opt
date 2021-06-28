#!/usr/bin/env python3

import numpy as np
import copy
import matscipy.calculators.eam.calculator
import shelve
import os
import sys
import shutil
from shutil import copyfile
import time


# Returns absolute value of distance from experimental value
def find_peak():
    try:
        os.system('./ns_analyse 32Cu.energies -M 0 -n 600 -D 5 > pf')
    except:
        print('Error calling ns_analyse')
        exit()
    T,Cp = numpy.loadtxt('pf', usecols=(0,3), unpack=True)
    index_PT =  numpy.where(Cp == numpy.amax(Cp))
    print(int(T[index_PT]))
    return(abs(experimental - int(T[index_PT])))

#Generates potential from list of parameters to modify original potential by. [embedded, density, potential]
def gen_poten(param):
    potential = matscipy.calculators.eam.io.read_eam('Cu_Zhou04.eam.alloy', 'eam/alloy')
    matscipy.calculators.eam.io.write_eam(potential[0], potential[1], (potential[2] * param[0]), (potential[3] * param[1]), (potential[4] *  param[2]), 'test.eam.alloy','eam/alloy')

# Starts nested sampling calculation
def call_ns_run(param):
    shutil.copytree('run', str(call_ns_run.counter))
    os.chdir(str(call_ns_run.counter))
    try:
        os.system('gen-ns') 
    except:
        print('Error Calling gen-ns')
        exit()
    gen_poten(param)
    print('Calling ns')
    runtime = time.time()
    try:
        os.system('srun ./ns_run < 32Cu.input')
    except:
        print('Error calling ns_run')
        exit()
    print('NS finished')
    runtime = (time.time() - runtime) / 3600
    h_time = open('times', 'a')
    h_time.write(str(runtime))
    h_time.close()

    score = find_peak()
    os.chdir('cd ../')
    call_ns_run.counter += 1
    return score 


call_ns_run.counter = 0

# NM parameters
step = 0.05
no_improv_thr = 20
n_improv_break= 3
max_iter = 15
alpha = 1.
gamma = 2.
rho = -0.5
sigma = 0.5

# Starting parameters
x_start = np.array([1.0, 1.0, 1.0])
experimental = 1400

# Progress file
h_progress = open('progress', 'a')
h_progress.write('NS count:  Stage:  Parameters:   Score:')



# init
dim = len(x_start)
prev_best = call_ns_run(x_start)
no_improv = 0
res = [[x_start, prev_best]]

for i in  range(dim):
    x = copy.copy(x_start)
    x[i] = x[i] + step
    gen_poten(x)
    score = call_ns_run(x)
    h_progress.write(str(call_ns_run.counter), 'Init', str(x), str(score))
    res.append([x, score])

# simplex iter
iters = 0
while 1:
    # order
    res.sort(key=lambda x: x[1])
    best = res[0][1]

    # break after max_iter
    if max_iter and iters >= max_iter:
        print('Reached max number of iterations')
        print( res[0])
        exit()
    iters += 1

    # break after no_improv_break iterations with no improvement
    print( '...best so far:', str(best))
    h_best = open('track_best', 'a')
    h_best.write('best so far:' + str(res[0]) + '\n')
    h_best.close()

    h_track = open('track', 'a')
    h_track.write(str(iters) + ': \n' + str(res))
    if best < prev_best - no_improve_thr:
        no_improv = 0
        prev_best = best
    else:
        no_improv += 1

    if no_improv >= no_improv_break:
        print('Reached max number of iterations without a significant improvement')
        print(res[0])
        exit()

    # centroid
    x0 = [0.] * dim
    for tup in res[:-1]:
        for i, c in enumerate(tup[0]):
            x0[i] += c / (len(res)-1)

    # reflection
    xr = x0 + alpha*(x0 - res[-1][0])
    gen_poten(xr)
    rscore = call_ns_run(xr)
    h_progress.write(str(call_ns_run.counter), 'Reflection', str(xr), str(score))
    if res[0][1] <= rscore < res[-2][1]:
        del res[-1]
        res.append([xr, rscore])
        continue

    # expansion
    if rscore < res[0][1]:
        xe = x0 + gamma*(x0 - res[-1][0])
        gen_poten(xe)
        escore = call_ns_run(xe)
        h_progress.write(str(call_ns_run.counter), 'Expansion', str(xe), str(score))
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
    gen_poten(xc)
    cscore = call_ns_run(xc)
    h_progress.write(str(call_ns_run.counter), 'Contraction', str(xc), str(score))
    if cscore < res[-1][1]:
        del res[-1]
        res.append([xc, cscore])
        continue

    # reduction
    x1 = res[0][0]
    nres = []

    #for tup in res:
    for tup in res: # might be a bug here from the loop type change, check
        redx = x1 + sigma*(tup[0] - x1)
        gen_poten(redx)
        score = call_ns_run(redx)
        h_progress.write(str(call_ns_run.counter), 'Reduction', str(redx), str(score))
        nres.append([redx, score])
    res = nres






