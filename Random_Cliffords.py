import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
import random

H = (np.array([[1, 1], [1, -1]])/np.sqrt(2)).astype(complex)
S = np.array([[1, 0], [0, 1j]]).astype(complex)
V = H @ S @ H @ S
W = H @ S

class Random_Clifford_Tableau: 
    '''
    Class for storing and evolving Clifford tableau for a qubit-chain
    via Aaronson-Gottesman algorithm
    '''
    def __init__(self, L = 10, D = 0, nT = 10):
        '''
        L: length of qubit chain
        D: number of single-qubit defect chains
        nT: number of Floquet periods
        '''
        assert (L > 0 and L < 1000)
        assert (D >= 0 and D <= 10)
        assert (nT >= 1 and nT <= 10)

        self.L = L
        self.D = D
        self.nT = nT
        
        # Integer stabiliser tableau for Pauli's and phases
        self.tableau = np.zeros((1, 2*self.L + 1)).astype(int)# Initial tableau with only identities
        # Weights of stabilisers of a superp
        self.weights = np.ones((1))
        self.tableau[:, -1] = np.ones(len(self.tableau[:, 0])).astype(int) # Set initial Pauli phase to 1

        self.z_indices = list(np.arange(0, self.L))
        self.x_indices = list(np.arange(self.L, 2*self.L))

    def Hadamard(self, i):
        '''
        Evolve tableau with Hadamard on qubit i
        '''
        
        # Swap Xi coloumn with Zi coloumn

        # [0, L-1]: z_indices
        # [L, 2L-1]: x_indices
        
        # Swap the ith z-coloumn with the ith x coloumnt
        self.tableau[:, [i, self.L+i]] = self.tableau[:, [self.L+i, i]]
        
        # Update phase of each coloumn
        self.tableau[:, 2*self.L] = self.tableau[:, 2*self.L]^(self.tableau[:, self.L+i]*self.tableau[:, i])
    
    def Phase(self, i):
        '''
        Evolve tableau with Phase gate on qubit i
        '''
        # [0, L-1]: z_indices
        # [L, 2L-1]: x_indices
        # Perform XOR on the ith Z coloumn and the ith X coloumn
        self.tableau[:, i] = self.tableau[:, i]^self.tableau[:, self.L+i]

        # Update phase of each coloumn
        self.tableau[:, 2*self.L] = self.tableau[:, 2*self.L]^(self.tableau[:, self.L+i]*self.tableau[:, i])
        
        
    def CNOT(self, i, j):
        '''
        Evolve tableau with CNOT with qubit i as control
        and qubit j as target
        '''

        self.tableau[:, i] = self.tableau[:, i]^self.tableau[:, j]
        self.tableau[:, self.L+j] = self.tableau[:, self.L+j]^self.tableau[:, self.L+i]
        # Update phase of each coloumn
        self.tableau[:, 2*self.L] = self.tableau[:, 2*self.L]^(self.tableau[:, self.L+i]*self.tableau[:, j])*\
                                    (self.tableau[:, self.L+j]^self.tableau[:, i]^np.ones(len(self.tableau[:, i])).astype(int))
    
    def X(self, i):
        '''
        Apply X-evolution on qubit i
        '''
        self.tableau[:, 2*self.L] = self.tableau[:, 2*self.L]^self.tableau[:, i]
    def Y(self, i):
        '''
        Apply iY-evolution on qubit i
        '''
        self.tableau[:, 2*self.L] = self.tableau[:, 2*self.L]^self.tableau[:, self.L+i]^self.tableau[:, i]
    def Z(self, i):
        '''
        Apply Z-evolution on qubit i
        '''
        self.tableau[:, 2*self.L] = self.tableau[:, 2*self.L]^self.tableau[:, self.L+i]
    def PH(self, i):
        '''
        Apply Phase+Hadamard on qubit i
        '''
        temp = np.copy(self.tableau[:, i])
        self.tableau[:, i] = self.tableau[:, i]^self.tableau[:, self.L+i]
        self.tableau[:, self.L+i] = temp

    def HP(self, i):
        '''
        Apply Hadamard+Phase on qubit i
        '''

        temp = np.copy(self.tableau[:, self.L+i])
        self.tableau[:, self.L+i] = self.tableau[:, self.L+i]^self.tableau[:, self.i]
        self.tableau[:, i] = temp

        self.tableau[:, 2*self.L] = self.tableau[:, 2*self.L]^self.tableau[:, self.L+i]

    def apply_random_clifford(self, i, j):
        '''
        Apply randomly drawn 2-qubit Clifford gate
        between sites i, j
        '''

        # Apply C1C2 random single-qubit Cliffords
        group1 = ["H", "I"]
        group2 = ["HP", "PH", "I"]
        group3 = ["X", "Y", "Z", "I"]

        rand = np.randint(0, 3)



        return 0
    
    def clifford_layer(self):
        '''
        Apply a layer of pairwise, random Clifford gates 
        '''
    def defect_layer(self):
        '''
        Apply a layer of randomly chosen defect gates
        '''
    def progress_time(self):
        '''
        Apply full-period of Floquet-evolution
        '''

    
