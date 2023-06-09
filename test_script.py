import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
import random

from Random_Cliffords import Random_Clifford_Tableau

tab = Random_Clifford_Tableau(L=200)
print(tab.tableau)

# Initialize with X on first qubit
tab.tableau[0, tab.L] = 1
print(tab.get_Pauli(0))
# Apply Brickwork evolution for 100 times
paulis = []
stats = []
T = 100
for t in range(0, T):
    tab.apply_brick_layer()
    pauli, stat = tab.get_Pauli(0, stats=True)
    print(pauli)
    print(stat)
    paulis.append(pauli)
    stats.append(stat)

print(stats)

stats= np.array(stats)
print(stats)
plt.plot(np.arange(0, T), stats[:,0]/tab.L, label='I')
plt.plot(np.arange(0, T), stats[:,1]/tab.L, label='X')
plt.plot(np.arange(0, T), stats[:,2]/tab.L, label='Y')
plt.plot(np.arange(0, T), stats[:,3]/tab.L, label='Z')
plt.xlabel('Periods')
plt.ylabel('Fraction of gates on chain')
plt.legend()
plt.savefig('stats.jpg')