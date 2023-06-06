import numpy as np
import scipy as sp
import matplotlib.pyplot as plt

H = (np.array([[1, 1], [1, -1]])/np.sqrt(2)).astype(complex)
S = np.array([[1, 0], [0, 1j]]).astype(complex)
V = H @ S @ H @ S
W = H @ S

class Random_Clifford_Tableau: 
    '''
    Class for storing and evolving Clifford tableau for a qubit-chain
    via Aaronson-Gottesman algorithm
    '''
    
    def __init__(L = 10, D = 0, nT = 10):
        '''
        L: length of qubit chain
        p: probability of single-qubit defect gate
        nT: number of Floquet periods
        '''
        assert L >= 0 and L < 1000
        assert D >= 0 and D <= 10
        assert nT >= 1 and nT < 10

        self.L = L
        self.D = D
        self.nT = nT

        self.tableau = np.zeros(2*L + 2) # Initial tableau with only identities
        self.tableau[-2] = 1 # Set initial Pauli phase to 1
        self.tableau[-1] = 1 # Set initial Pauli weight to 1

        return 
    
    def apply_random_clifford(i, j):
        '''
        Apply randomly drawn 2-qubit Clifford gate
        between sites i, j
        '''
        return
    
    def clifford_layer():
        '''
        Apply a layer of pairwise, random Clifford gates 
        '''
    def defect_layer():
        '''
        Apply a layer of randomly chosen defect gates
        '''
    def progress_time():
        '''
        Apply full-period of Floquet-evolution
        '''

    
