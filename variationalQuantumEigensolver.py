from qiskit import *
from qiskit_aer import *

def cvarExpectedValue(counts, shots, alpha):
    '''
    Calculates the Conditional Value at Risk of a Hamiltonian given the counts from a quantum circuit
    
    Parameters:
    counts (dict): A dictionary containing the counts from a quantum circuit
    shots (int): The number of shots
    alpha (float): The confidence level

    Returns:
    value (float): The Conditional Value at Risk of the Hamiltonian
    '''
    values = []
    
    for bitstring, count in counts.items():
        sign = (-1) ** bitstring.count('1')
        values.extend([sign] * count)
    
    values.sort()

    cutoffIndex = int(alpha * shots)

    cvarSum = sum(values[:cutoffIndex])

    cvarValue = cvarSum / cutoffIndex

    return cvarValue

def ParametricCircuitVQE(circ, n, nLayers, theta):
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
    for layer in range(nLayers):
        for qubit in range(n):
            circ.ry(theta[counter], qubit)
            counter += 1
            if qubit < n - 1:
                circ.cx(qubit, qubit + 1)
    return circ

def NParaVQE(n):
    '''
    Calculates the number of parameters needed for the VQE algorithm
    
    Parameters:
    n (int): The number of qubits
    
    Returns:
    n * k (int): The number of parameters
    k (int): The number of layers
    '''
    k = 1
    while n * k <= n ** 2: 
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

def cvarVQE(theta, quantumParameters):
    '''
    Calculates the Conditional Value at Risk of the Hamiltonian given the parameters theta
    
    Parameters:
    theta (list): A list of parameters
    quantumParameters (dict): A dictionary containing quantum parameters
    
    Returns:
    ev (float): The Conditional Value at Risk of the Hamiltonian
    '''
    n = quantumParameters['nVQE']
    nLayers = quantumParameters['nLayersVQE']
    h = quantumParameters['hVQE']
    shots = quantumParameters['shots']
    alpha = quantumParameters['cvarAlpha']
    ev = 0
    for key in h:
        if key != '0' * n:
            circ = QuantumCircuit(n, n)
            circ = ParametricCircuitVQE(circ, n, nLayers, theta)
            for j in range(n):
                if key[j] == '3':
                    circ.measure(j, j)
            counts = Sim(circ, shots)
            value = cvarExpectedValue(counts, shots, alpha)
            ev += h[key] * value
        else:
            ev += h[key]
    return ev

def thetaToSolutionVQE(theta, parameters):
    '''
    Converts the parameters theta to a solution for the VQE algorithm
    
    Parameters:
    theta (list): A list of parameters
    parameters (dict): A dictionary containing wind farm parameters
    
    Returns:
    solution (list): A list representing the wind farm layout solution
    '''
    circ = QuantumCircuit(parameters['nVQE'], parameters['nVQE'])
    circ = ParametricCircuitVQE(circ, parameters['nVQE'], parameters['nLayersVQE'], theta)
    for j in range(parameters['nVQE']):
        circ.x(j)
        circ.measure(j, j)
    counts = Sim(circ, parameters['shots'])
    
    bestKey = max(counts, key=counts.get)
    return [int(i) for i in bestKey]
