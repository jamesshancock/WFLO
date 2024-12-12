from qiskit import *
from qiskit_aer import *
import numpy as np
from extras import *

def ExpectedValue(counts, shots):
    '''
    Calculates the expected value of a Hamiltonian given the counts from a quantum circuit
    
    Parameters:
    counts (dict): A dictionary containing the counts from a quantum circuit
    shots (int): The number of shots
    
    Returns:
    value (float): The expected value of the Hamiltonian
    '''
    n = len(list(counts.keys())[0])
    value = 0
    for key in counts:
        sign = 1
        for i in key:
            if i == '1':
                sign *= -1
        value += sign * counts[key] / shots
    return value

def ParametricCircuitPCE(circ, n, nLayers, theta):
    '''
    Generates the parametric circuit for the VQE algorithm
    
    Parameters:
    circ (QuantumCircuit): The quantum circuit
    n (int): The number of qubits
    nLayers (int): The number of layers
    theta (list): A list of parameters
    
    Returns:
    circ (QuantumCircuit): The quantum circuit with the param
    '''
    counter = 0
    for layer in range(int(nLayers/2)):
        for qubit in range(n):
            circ.ry(theta[counter], qubit)
            circ.ry(theta[counter + 1], qubit)           
            if qubit < n - 1:
                circ.cx(qubit, qubit + 1)
            counter += 2
    return circ

def NParaPCE(n):
    '''
    Calculates the number of parameters needed for the VQE algorithm
    
    Parameters:
    n (int): The number of qubits
    
    Returns:
    n * k (int): The number of parameters
    k (int): The number of layers
    '''
    k = 1
    while n * k <= 2**(n+1) + 4*n: 
        k += 1
    k -= 1
    return n * k, k

def Sim(circ, shots):
    '''
    Simulates the quantum circuit - returns the counts
    
    Parameters:
    circ (QuantumCircuit): The quantum circuit
    shots (int): The number of shots
    
    Returns:
    counts (dict): A dictionary containing
    '''
    backendSim = Aer.get_backend('qasm_simulator')
    jobSim = backendSim.run(transpile(circ, backendSim), shots=shots)
    resultSim = jobSim.result()
    return resultSim.get_counts(circ)

def ExpectedValueForKey(theta, quantumParameters, key):
    '''
    Calculates the expected value of the Hamiltonian given the parameters theta
    
    Parameters:
    theta (list): A list of parameters
    quantumParameters (dict): A dictionary containing quantum parameters
    
    Returns:
    ev (float): The expected value of the Hamiltonian
    '''
    n = quantumParameters['nPCE']
    nLayers = quantumParameters['nLayersPCE']
    shots = quantumParameters['shots']
    circ = QuantumCircuit(n, n)
    circ = ParametricCircuitPCE(circ, n, nLayers, theta)
    for j in range(n):
        if key[j] == '3':
            circ.measure(j, j)
        elif key[j] == '1':
            circ.h(j)
            circ.measure(j, j)
        elif key[j] == '2':
            circ.sdg(j)
            circ.h(j)
            circ.measure(j, j)
    counts = Sim(circ, shots)
    return np.tanh(quantumParameters['talpha']*ExpectedValue(counts, shots))

def PCE(theta, parameters, keys, Jprime, hPCE):
    '''
    Calculates the PCE cost function
    
    Parameters:
    theta (list): A list of parameters
    parameters (dict): A dictionary containing parameters
    keys (list): A list of keys
    Jprime (numpy.ndarray): The upper triangular matrix
    
    Returns:
    cost (float): The cost of the PCE
    '''
    store = []
    cost = 0
    for key in keys:
        store.append(ExpectedValueForKey(theta, parameters, key))
    for i in range(Jprime.shape[0]):
        for j in range(i+1,Jprime.shape[0]):
            cost += Jprime[i,j]*store[i]*store[j]
        cost += hPCE[i]*store[i] 
    return cost

def thetaToSolutionPCE(theta, parameters, keys):
    '''
    Converts the parameters theta to a solution for the PCE algorithm
    
    Parameters:
    theta (list): A list of parameters
    parameters (dict): A dictionary containing wind farm parameters
    keys (list): A list of keys
    
    Returns:
    solution (list): A list representing the wind farm layout solution
    '''
    store = []
    for key in keys:
        EV = ExpectedValueForKey(theta, parameters, key)
        if EV > 0:
            store.append(1)
        else:
            store.append(0)    
    return store
