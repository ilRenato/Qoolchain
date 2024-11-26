from libcpp.vector cimport vector
from libcpp.utility cimport pair
from libcpp.map cimport map
from libcpp cimport bool
from libcpp.list cimport list as cpplist
from libcpp.string cimport string

cdef extern from "PreQubo.h":
    pair[vector[pair[vector[vector[int]], cpplist[int]]], map[int, bool]] PreQubo(vector[vector[int]]& f, int &offset, string &filename)
    int getOffset()
    float getLowerBound()
    vector[pair[vector[pair[vector[vector[int]], cpplist[int]]], map[int, bool]]] ShannonDecomposition(vector[vector[int]] &Q, int &n_iterations, int &flag)

def PreQuboWrap(Q, offset, filename):
    decompositions, persistencies = PreQubo(Q, offset, filename)
    lower_bound = getLowerBound()
    subqubos = []
    varnames = []
    for q, v in decompositions:
        subqubos.append(q)
        varnames.append(v)
    return subqubos, varnames, persistencies, lower_bound

def ShannonDecompositionWrap(Q, n_iterations, flag):
    qubos_list = ShannonDecomposition(Q, n_iterations, flag)
    return qubos_list
