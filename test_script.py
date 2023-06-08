import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
import random

from Random_Cliffords import Random_Clifford_Tableau

tab = Random_Clifford_Tableau(L=5)
print(tab.tableau)

tab.tableau[0, 0] = 1
print(tab.tableau)
tab.Hadamard(0)
print(tab.tableau)
tab.Phase(0)
print(tab.tableau)
tab.CNOT(0, 1)
print(tab.tableau)