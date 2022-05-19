#!/usr/bin/env python3

import copy
import os
import shutil
import time
import matscipy.calculators.eam.calculator
import subprocess
from config import *


# Returns absolute value of distance from experimental value

def find_peak(prefix, target):
    if os.path.exists(str(prefix) + '_pt_temp'):
        pt_temp = int(np.loadtxt(str(prefix) + '_pt_temp'))
        return abs(target - (pt_temp + finite_size_offset))

    else:
        subprocess.run(str(NS_location) + 'ns_analyse ' + str(prefix) + '.energies -M 1 -n 1000 -D 5 > ' + str(prefix)
                       + '.dat', check=True, shell=True)

        t, cp = np.loadtxt((str(prefix) + '.dat'), usecols=(0, 4), unpack=True)
        index_of_phase_transition = np.where(cp == np.amax(cp))
        # If multiple temperatures have the same value for heat capacity, will select the lowest temperature of those
        # with the same heat capacity
        index_of_phase_transition = index_of_phase_transition[0]
        pt_temp = int(t[index_of_phase_transition])

        h_pt_temp = open(str(prefix) + '_pt_temp', 'a')
        h_pt_temp.write(str(pt_temp))
        h_pt_temp.close()

        return abs(target - (pt_temp + finite_size_offset))


# Function to multiply values in the potential with a linear line across the radius.
# Create line with the same number of value as the potential and multiply against the potential
# If statement to deal with if parameters are the same resulting in a flat line - Line spacing calculated using
# difference between parameters so need another method when they're the same.

def gen_poten(param):
    potential = matscipy.calculators.eam.io.read_eam(potential_name, potential_type)
    if method == 'line':

        graph_1 = potential[2]
        graph_2 = potential[3]
        graph_3 = potential[4]

        if param[1] == param[0]:
            line_1 = np.full(np.size(graph_1), param[0])
        else:
            line_1 = np.arange(param[0], param[1], ((param[1] - param[0]) / np.size(graph_1)))

        if param[2] == param[3]:
            line_2 = np.full(np.size(graph_2), param[2])
        else:
            line_2 = np.arange(param[2], param[3], ((param[3] - param[2]) / np.size(graph_2)))

        if param[4] == param[5]:
            line_3 = np.full(np.size(graph_3), param[4])
        else:
            line_3 = np.arange(param[4], param[5], ((param[5] - param[4]) / np.size(graph_3)))

        new_graph_1 = (np.multiply(graph_1, line_1))
        new_graph_2 = (np.multiply(graph_2, line_2))
        new_graph_3 = (np.multiply(graph_3, line_3))

        matscipy.calculators.eam.io.write_eam(potential[0], potential[1], new_graph_1, new_graph_2, new_graph_3,
                                              'test.eam.alloy', 'eam/alloy')
    elif method == 'factor':

        matscipy.calculators.eam.io.write_eam(potential[0], potential[1], (potential[2] * param[0]),
                                              (potential[3] * param[1]), (potential[4] * param[2]), 'test.eam.alloy',
                                              'eam/alloy')
    else:
        print('No method defined')
        exit()


# Starts nested sampling calculation

def call_ns_run(param):
    # See if folder exists already, if it does we've restarted, and we want to restart the ns calculation
    global call_ns_run_counter
    iter_score = 0
    if not os.path.isdir(str(call_ns_run_counter)):
        shutil.copytree('run', str(call_ns_run_counter))
    os.chdir(str(call_ns_run_counter))
    gen_poten(param)
    print('Calling ns')
    runtime = time.time()

    for j in range(len(prefixes)):
        subprocess.run('srun ' + str(NS_location) + 'ns_run < ' + str(prefixes[j]) + '.input', check=True, shell=True)
        iter_score += iter_score + find_peak(prefixes[j], targets[j])

    print('NS finished')
    runtime = (time.time() - runtime) / 3600
    h_time = open('times', 'a')
    h_time.write(str(runtime))
    h_time.close()

    os.chdir('../')
    call_ns_run_counter += 1
    return iter_score


def write_progress_file(stage=None, current_parameter=None, current_score=None):
    if not os.path.exists('progress'):
        h_progress = open('progress', 'a')
        h_progress.write('NS count:  Stage:  Parameters:   Score: \n')
        h_progress.close()
    else:
        h_progress = open('progress', 'a')
        h_progress.write(str(call_ns_run_counter) + ' ' + str(stage) + str(current_parameter) + ' ' + str(current_score)
                         + '\n')
        h_progress.close()


call_ns_run_counter = 0
write_progress_file()

# initialization
dim = len(x_start)
prev_best = call_ns_run(x_start)
write_progress_file('Init', x_start, prev_best)
no_improv = 0
res = [[x_start, prev_best]]
iters = 0

for i in range(dim):
    x = copy.copy(x_start)
    x[i] = x[i] + step
    score = call_ns_run(x)
    write_progress_file('Init', x, score)
    res.append([x, score])

# simplex iter
while 1:
    # order
    res.sort(key=lambda o: o[1])
    best = res[0][1]

    # break after max_iter
    if max_iter and iters >= max_iter:
        print('Reached max number of iterations')
    iters += 1

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
    xr = x0 + alpha * (x0 - res[-1][0])
    rscore = call_ns_run(xr)
    write_progress_file('Reflection', xr, rscore)
    if res[0][1] <= rscore < res[-2][1]:
        del res[-1]
        res.append([xr, rscore])
        continue

    # expansion
    if rscore < res[0][1]:
        xe = x0 + gamma * (x0 - res[-1][0])
        escore = call_ns_run(xe)
        write_progress_file('Expansion', xe, escore)
        if escore < rscore:
            del res[-1]
            res.append([xe, escore])
            continue
        else:
            del res[-1]
            res.append([xr, rscore])
            continue

    # contraction
    xc = x0 + rho * (x0 - res[-1][0])
    cscore = call_ns_run(xc)
    write_progress_file('Contraction', xc, cscore)
    if cscore < res[-1][1]:
        del res[-1]
        res.append([xc, cscore])
        continue

    # reduction
    x1 = res[0][0]
    nres = []

    # for tup in res:
    for tup in res:  # might be a bug here from the loop type change, check
        redx = x1 + sigma * (tup[0] - x1)
        score = call_ns_run(redx)
        write_progress_file('Reduction', redx, score)
        nres.append([redx, score])
    res = nres
