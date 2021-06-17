#!/usr/bin/env python3

import numpy as np
import copy
import matscipy.calculators.eam.calculator

poten_intro =[]


def split_poten(poten):
    values=[]
    lines =[]

    n=0
    with open(str(poten)) as file:
        func_interatomic = [] 
        func_embedding   = []
        func_density     = []
        for line in file:
            if n < 6 :
                n += 1
                poten_intro.append(line)
                continue
            line=line.split()
            for c in range(len(line)):
                values.append(line[c])
    for n in range(len(values)):
            if n < 2000:
                func_embedding.append(values[n])
            elif n< 4000:
                func_interatomic.append(values[n])
            else :
                func_density.append(values[n])
    return(func_embedding, func_interatomic, func_density)


def write_poten(vertex, filename):
    h_poten = open(str(filename), 'w')
   
    for x in range(len(poten_intro)):
        h_poten.write(poten_intro[x])

    for x in range(len(embedding[vertex])):
        if (x  % 5) == 0 and x != 0 :
            h_poten.write('\n')
        h_poten.write(str(embedding[vertex]) + ' ')

    for x in range(len(interatomic[vertex])):
        if (x  % 5) == 0 :
            h_poten.write('\n')
        h_poten.write(str(interatomic[vertex]) + ' ')

    for x in range(len(density[vertex])):
        if (x  % 5) == 0 :
            h_poten.write('\n')
        h_poten.write(str(density[vertex]) + ' ')
    h_poten.close()

def find_peak():
    os.system('gen-ns')
    os.system('./ns_analyse 64Cu.energies -M 0 -n 600 -D 5 > pf')
    T,Cp = numpy.loadtxt('pf', usecols=(0,3), unpack=True)
    index_PT =  numpy.where(Cp == numpy.amax(Cp))
    print(int(T[index_PT]))
    return(int(T[index_PT]))

with open('Cu_Zhou04.eam.alloy') as file:
    n = 0
    for line in file:
        if n < 6 :
            n += 1
            poten_intro.append(line)


orig_embedding, orig_interatomic, orig_density = split_poten('Cu_Zhou04.eam.alloy')

#embedding = np.array([orig_embedding]*4)
interatomic = np.array([orig_interatomic]*4)
density = np.array([orig_density]*4)


potential = matscipy.calculators.eam.io.read_eam('Cu_Zhou04.eam.alloy', 'eam/alloy')

embedding = np.array([potential[2],]*4)
interatomic = np.array([potential[4],]*4)
density =np.array([potential[3],]*4) 
