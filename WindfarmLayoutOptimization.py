# Optimizers: COByLA, Bayesian (NEW kernel !!)
'''
Packages
'''
import numpy as np
import matplotlib.pyplot as plt

'''
Modules
'''
from windfarmQUBO import *
from variationalQuantumEigensolver import *
from pauliCorrelationEncoding import *
from solvers import *
from extras import *
#[[10*k, 12, 1/36] for k in range(36)]

parameters = {
           'len_grid': 4, 
           'D': [[0,12,1]], 
           'x': 1, 
           'r': 1.5, 
           'm': 2,    # Number of turbines
           'E': 1.0,    # Proximity threshold
           'lam1': 1e4, # Number constraint
           'lam2': 3e3, # Proximity constraint
           'shots': 3000,
           'cvarAlpha': 0.8,
           'solver': 'PCE',
           'method': 'COBYLA',
           'nSamples': 3,
           'talpha': 10,
           'L': 15,
           'tol': 1,
           'bayesSigma': 1000.0,
           'bayesGamma': 0.01,
           'bayesIters': 100,
           'ExhaustiveCheck': False,
           }
# TEST ZONE
# =============================================================================


# ============================================================================= #
if parameters['ExhaustiveCheck']:
    optEnergy = printFaff(parameters)
    optSpinCost, optSpinVec = spinExhaustiveCheck(WindfarmQ(parameters), parameters)
    print("The optimal spin energy is", optSpinCost)
    print("One optimal solution is", optSpinVec)


print("----------------------------")
print("Quantum Solvers")
print("----------------------------")
Energies = []
labels = []


print("Solver:", parameters['solver'])
print("Method:", parameters['method'])

solution, history, timeTaken = quantumSolver(parameters)

print("Solution:")
print(solution)
print("Grid:")
print(solutionToGrid(solution, parameters))












