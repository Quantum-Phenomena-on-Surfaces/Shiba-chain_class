# -*- coding: utf-8 -*-
"""
Created on Mon Feb  1 12:12:22 2021

@author: crisl
"""

import numpy as np
from numpy.linalg import inv
from scipy.sparse import identity
from scipy.sparse import csc_matrix
#from scipy.sparse.linalg import inv
from scipy.sparse.linalg import spsolve

#solve Dyson's equation
def Dyson_eq(Go , Self , N_x, N_y, layers):    
    
    Id = np.identity(4 * N_y * N_x * layers)
    #Id = identity(4 * N_y * N_x * layers)
    #Id = csc_matrix(Id1)
     
    S = Id - np.dot(Go, Self)
    #matrx_inv = spsolve(S, Id)
    matrx_inv = inv(Id - np.dot(Go, Self))##### + o -?????
    gg = np.dot(matrx_inv , Go)      
    
    #covert back to array
    #gg = gg.toarray()
    
    return(gg)
