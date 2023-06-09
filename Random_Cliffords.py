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
        self.tableau[:, self.L+i] = self.tableau[:, self.L+i]^self.tableau[:, i]
        self.tableau[:, i] = temp

        self.tableau[:, 2*self.L] = self.tableau[:, 2*self.L]^self.tableau[:, self.L+i]
    def apply_random_Pauli(self, i):
        
        rand = random.randint(0, 4)
        if (rand == 0):
            return
        elif (rand == 1):
            self.X(i)
        elif (rand == 2):
            self.Y(i)
        elif (rand == 3):
            self.Z(i)

    def apply_random_H(self, i):

        rand = random.randint(0, 2)
        if (rand == 0):
            return
        elif (rand == 1):
            self.Hadamard(i)
    
    def apply_random_V(self, i):

        rand = random.randint(0, 3)
        if (rand == 0):
            return
        elif (rand == 1):
            self.HP(i)
        elif (rand == 2):
            self.PH(i)

    def apply_random_clifford1(self, i):
        '''
        Apply randomly drawn 1-qubit Clifford gate on site i
        '''
        self.apply_random_H(i)
        self.apply_random_V(i)
        self.apply_random_Pauli(i)
    
    def apply_random_clifford2(self, i, j):
        '''
        Apply randomly drawn 2-qubit Clifford on site i and j
        '''

        self.apply_random_clifford1(j)
        self.apply_random_clifford1(i)

        rand = random.randint(0, 20)
        if (rand == 0):
            return
        elif (rand >=1 and rand <= 9):
            self.CNOT(i, j)
            self.apply_random_V(j)
            self.apply_random_V(i)
        elif (rand >= 10 and rand <= 18):
            self.CNOT(i, j)
            self.CNOT(j, i)
            self.apply_random_V(j)
            self.apply_random_V(i)
        elif (rand == 19):
            self.CNOT(i, j)
            self.CNOT(j, i)
            self.CNOT(i, j)

    def get_Pauli(self, i, stats=False):
        '''
        Calculate the Pauli string for each row
        '''
        # ab with a: Z component, b X component
        # 00 --> I
        # 01 --> X
        # 10 --> Z
        # 11 --> Y

        mapping = [['I', 'X'], ['Z', 'Y']]
        gate_counts = [[0, 0], [0, 0]]
        p_string = []

        for idx in range(0, self.L):
            p_string.append(mapping[self.tableau[i, idx]][self.tableau[i, self.L+idx]])
            gate_counts[self.tableau[i, idx]][self.tableau[i, self.L+idx]] += 1
        if (stats):
            nI = gate_counts[0][0]
            nX = gate_counts[0][1]
            nY = gate_counts[1][1]
            nZ = gate_counts[1][0]

            return p_string, [nI, nX, nY, nZ]

        return p_string


    def apply_brick_layer(self):
        '''
        Apply nearest-neighbour random Clifford gates for whole chain
        '''
        for qubit_idx in range(0, self.L-1):
            self.apply_random_clifford2(qubit_idx, qubit_idx+1)
        #if (self.L%2 == 1):
        #   self.apply_random_clifford2(self.L-1, 0)
        for qubit_idx in range(1, self.L-1):
            self.apply_random_clifford2(qubit_idx, qubit_idx+1)
        #self.apply_random_clifford2(self.L-1, 0)
        #if (self.L%2 == 1):
        #    self.apply_random_clifford2(0, 1)


class Random_Clifford_Evolution(Random_Clifford_Tableau):
    def __init__(self, L, D, T):
        Random_Clifford_Tableau.__init__(self, L, D, T)


