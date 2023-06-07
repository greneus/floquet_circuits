import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
import random

class Random_Clifford_Tableau: 
    '''
    Class for storing and evolving Clifford tableau for a qubit-chain
    via Aaronson-Gottesman algorithm
    '''
    
    def __init__(L = 10, D = 0, nT = 10):
        '''
        L: length of qubit chain
        D: number of single-qubit defect chains
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

        self.z_indices = list(np.arange(0, L))
        self.x_indices = list(np.arange(L, 2*L))

        return 

    def Hadamard(i):
        '''
        Evolve tableau with Hadamard on qubit i
        '''
        
        # Swap Xi coloumn with Zi coloumn

        self.tableau[:, [i, L+i]] = self.tableau[:, [L+i, i]]
        
        # Update phase of each coloumn
        self.tableau[:, 2L] ^= self.tableau[:, i]*self.tableau[:, L+i]
    
    def Phase(i):
        '''
        Evolve tableau with Phase gate on qubit i
        '''
        # [0, L-1]: z_indices
        # [L, 2L-1]: x_indices
        self.tableau[:, i] = np.bitwise_xor(self.tableau[:, i], self.tableau[:, L+i])

         # Update phase of each coloumn
        self.tableau[:, 2L] ^= self.tableau[:, i]*self.tableau[:, L+i]

    def CNOT(i, j):
        '''
        Evolve tableau with CNOT with qubit i as control
        and qubit j as target
        '''

        self.tableau[:, i] = np.bitwise_xor(self.tableau[:, i], self.tableau[:, j])
        self.tableau[:, L+j] = np.bitwise_xor(self.tableau[:, L+j], self.tableau[:, L+i])
        self.tableau[:, 2L] ^= self.tableau[:, i]*self.tableau[:, L+j]\
                               *np.bitwise_xor(\
                                    np.bitwise_xor(self.tableau[:, L+j], self.tableau[:, i]), \
                                    np.ones(len(self.tableau[:, i])))
        
    
    def apply_random_clifford(i, j):
        '''
        Apply randomly drawn 2-qubit Clifford gate
        between sites i, j
        '''

        # Apply C1C2 random single-qubit Cliffords


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

    
