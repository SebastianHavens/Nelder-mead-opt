#!/usr/bin/env python3

import numpy as np
import copy
import matscipy.calculators.eam.calculator

poten_intro =[]

def split_poten(poten):
    values=[]
    lines =[]
    embedding = []
    interatomic = []
    density = []

    n=0
    with open(str(poten)) as file:
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
                embedding.append(values[n])
            elif n< 4000:
                interatomic.append(values[n])
            else :
                density.append(values[n])
    return(embedding, interatomic, density)


def write_poten(vertex, filename):
    h_poten = open(str(filename), 'w')
   
    for x in range(len(poten_intro)):
        h_poten.write(poten_intro[x])

    for x in range(len(vertex.embedding)):
        if (x  % 5) == 0 and x != 0 :
            h_poten.write('\n')
        h_poten.write(str(vertex.embedding[x]) + ' ')

    for x in range(len(vertex.interatomic)):
        if (x  % 5) == 0 :
            h_poten.write('\n')
        h_poten.write(str(vertex.interatomic[x]) + ' ')

    for x in range(len(vertex.density)):
        if (x  % 5) == 0 :
            h_poten.write('\n')
        h_poten.write(str(vertex.density[x]) + ' ')
    h_poten.close()

def find_peak():
    os.system('gen-ns')
    os.system('./ns_analyse 64Cu.energies -M 0 -n 600 -D 5 > pf')
    T,Cp = numpy.loadtxt('pf', usecols=(0,3), unpack=True)
    index_PT =  numpy.where(Cp == numpy.amax(Cp))
    print(int(T[index_PT]))
    return(int(T[index_PT]))

class vertex :
    T_PT = 9999
    def __init__(vertex, poten):
        vertex.embedding, vertex.interatomic, vertex.density = split_poten(str(poten))




vertex_original = vertex('Cu_Zhou04.eam.alloy')


vertex_0 = copy.deepcopy(vertex_original)
vertex_1 = copy.deepcopy(vertex_original)
vertex_2 = copy.deepcopy(vertex_original)
vertex_3 = copy.deepcopy(vertex_original)

vertex_1.embedding = [float (i) * 1.05 for i in vertex_1.embedding]
vertex_2.interatomic = [float (i) * 1.05 for i in vertex_2.interatomic]
vertex_3.density = [float (i) * 1.05 for i in vertex_3.density]


#write_poten(vertex_original, 'foo.eam.alloy')
write_poten(vertex_0, 'vertex_0.eam.alloy')
write_poten(vertex_1, 'vertex_1.eam.alloy')
write_poten(vertex_2, 'vertex_2.eam.alloy')
write_poten(vertex_3, 'vertex_3.eam.alloy')


print(vertex.embedding)
#make a class which has a list of the three function values 
# easier to track that way
