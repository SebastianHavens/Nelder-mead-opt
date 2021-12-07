#!/usr/bin/env python3

import numpy as np
import copy
import matscipy.calculators.eam.calculator
import os
import sys
import shutil
import time
import shelve


# Returns absolute value of distance from experimental value
def find_peak():
    try:
        os.system('./ns_analyse 32Cu.energies -M 1 -n 1000 -D 5 > pf')
    except:
        print('Error calling ns_analyse')
        exit()
    T, Cp = np.loadtxt('pf', usecols=(0, 4), unpack=True)
    index_PT = np.where(Cp == np.amax(Cp))
    print(int(T[index_PT]))
    
    return (abs(experimental - (int(T[index_PT]) - 224  ))) #Calculate error from experimental value, we take away an amount dependent on the system size effects from out calciulated phase transition temperature..


# Generates potential from list of parameters to modify original potential by. [embedded, density, potential]
def gen_poten(param):
    potential = matscipy.calculators.eam.io.read_eam('Cu_Zhou04.eam.alloy', 'eam/alloy')
    matscipy.calculators.eam.io.write_eam(potential[0], potential[1], (potential[2] * param[0]),
                                          (potential[3] * param[1]), (potential[4] * param[2]), 'test.eam.alloy',
                                          'eam/alloy')
    print('Current parameters:')
    print(param)


# Starts nested sampling calculation
def call_ns_run(param):
    # See if folder exists already, if it does we've restarted and we want to restart the ns calculation
    if not os.path.isdir(str(call_ns_run.counter)):
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
    os.chdir('../')
    call_ns_run.counter += 1
    return score


def load_variables():
    cwd = os.getcwd()
    shelf = shelve.open(str(cwd + '/variables'))
    for key in shelf:
        globals()[key] = shelf[key]
    shelf.close


def save_variables():
    cwd = os.getcwd()
    filename = str(cwd) + '/variables'
    shelf = shelve.open(filename, 'n')
    for key in globals():
        try:
            shelf[key] = globals()[key]
        except:
            print('ERROR shleving: {0}'.format(key))
    shelf.close()


opts = [opt for opt in sys.argv[1:] if opt.startswith("-")]
if "-r" in opts:
    print('Restarting NM calculation')
    restart = True
else:
    restart = False

# NM parameters
step = 0.05
no_improve_thr = 20
no_improv_break = 3
max_iter = 15
alpha = 1.
gamma = 2.
rho = -0.5
sigma = 0.5

# Starting parameters
x_start = np.array([1.0, 1.0, 1.0])
experimental = 1400

stage = None
if restart == True:
    load_variables()
    # reflection
    if stage == 0:
        xr = x0 + alpha * (x0 - res[-1][0])
        rscore = call_ns_run(xr)
        h_progress = open('progress', 'a')
        h_progress.write(str(call_ns_run.counter) + ' Reflection ' + str(xr) + ' ' + str(rscore) + '\n')
        h_progress.close()
        if res[0][1] <= rscore < res[-2][1]:
            del res[-1]
            res.append([xr, rscore])
        else:
            stage = 1

    # expansion
    if stage == 1:
        if rscore < res[0][1]:
            xe = x0 + gamma * (x0 - res[-1][0])
            escore = call_ns_run(xe)
            h_progress = open('progress', 'a')
            h_progress.write(str(call_ns_run.counter) + ' Expansion ' + str(xe) + ' ' + str(escore) + '\n')
            h_progress.close()
            if escore < rscore:
                del res[-1]
                res.append([xe, escore])
            else:
                del res[-1]
                res.append([xr, rscore])
        else:
            stage = 2

    # contraction
    if stage == 2:
        xc = x0 + rho * (x0 - res[-1][0])
        cscore = call_ns_run(xc)
        h_progress = open('progress', 'a')
        h_progress.write(str(call_ns_run.counter) + ' Contraction ' + str(xc) + ' ' + str(cscore) + '\n')
        h_progress.close()
        if cscore < res[-1][1]:
            del res[-1]
            res.append([xc, cscore])
        else:
            stage = 3

    # reduction
    if stage == 3:
        x1 = res[0][0]
        nres = []

        # for tup in res:
        for tup in res:  # might be a bug here from the loop type change, check
            redx = x1 + sigma * (tup[0] - x1)
            score = call_ns_run(redx)
            h_progress = open('progress', 'a')
            h_progress.write(str(call_ns_run.counter) + ' Reduction ' + str(redx) + ' ' + str(score) + '\n')
            h_progress.close()
            nres.append([redx, score])
        res = nres

# Progress file
if restart == False:
    call_ns_run.counter = 0
    h_progress = open('progress', 'a')
    h_progress.write('NS count:  Stage:  Parameters:   Score:')
    h_progress.close()
    # init
    dim = len(x_start)
    prev_best = call_ns_run(x_start)
    no_improv = 0
    res = [[x_start, prev_best]]
    iters = 0

    for i in range(dim):
        x = copy.copy(x_start)
        x[i] = x[i] + step
        score = call_ns_run(x)
        h_progress = open('progress', 'a')
        h_progress.write(str(call_ns_run.counter) + ' Init ' + str(x) + ' ' + str(score) + '\n')
        h_progress.close()
        res.append([x, score])

# simplex iter
while 1:
    # order
    res.sort(key=lambda x: x[1])
    best = res[0][1]

    # break after max_iter
    if max_iter and iters >= max_iter:
        print('Reached max number of iterations')
        print(res[0])
        exit()
    iters += 1

    # break after no_improv_break iterations with no improvement
    print('...best so far:', str(best))
    h_best = open('track_best', 'a')
    h_best.write('best so far: ' + str(res[0]) + '\n')
    h_best.close()

    h_track = open('track', 'a')
    h_track.write(str(iters) + ': ' + str(res) + '\n')
    h_track.close()
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
            x0[i] += c / (len(res) - 1)

    # reflection
    stage = 0
    save_variables()
    xr = x0 + alpha * (x0 - res[-1][0])
    rscore = call_ns_run(xr)
    h_progress = open('progress', 'a')
    h_progress.write(str(call_ns_run.counter) + ' Reflection ' + str(xr) + ' ' + str(rscore) + '\n')
    h_progress.close()
    if res[0][1] <= rscore < res[-2][1]:
        del res[-1]
        res.append([xr, rscore])
        continue

    # expansion
    if rscore < res[0][1]:
        stage = 1
        save_variables()
        xe = x0 + gamma * (x0 - res[-1][0])
        escore = call_ns_run(xe)
        h_progress = open('progress', 'a')
        h_progress.write(str(call_ns_run.counter) + ' Expansion ' + str(xe) + ' ' + str(escore) + '\n')
        h_progress.close()
        if escore < rscore:
            del res[-1]
            res.append([xe, escore])
            continue
        else:
            del res[-1]
            res.append([xr, rscore])
            continue

    # contraction
    stage = 2
    save_variables()
    xc = x0 + rho * (x0 - res[-1][0])
    cscore = call_ns_run(xc)
    h_progress = open('progress', 'a')
    h_progress.write(str(call_ns_run.counter) + ' Contraction  ' + str(xc) + ' ' + str(cscore) + '\n')
    h_progress.close()
    if cscore < res[-1][1]:
        del res[-1]
        res.append([xc, cscore])
        continue

    # reduction
    stage = 3
    save_variables()
    x1 = res[0][0]
    nres = []

    # for tup in res:
    for tup in res:  # might be a bug here from the loop type change, check
        redx = x1 + sigma * (tup[0] - x1)
        score = call_ns_run(redx)
        h_progress = open('progress', 'a')
        h_progress.write(str(call_ns_run.counter) + ' Reduction  ' + str(redx) + ' ' + str(score) + '\n')
        h_progress.close()
        nres.append([redx, score])
    res = nres
