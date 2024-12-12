import math
import numpy as np

def Labelling(lenGrid):
    '''
    Creates a dictionary with the labels of the grid points as keys and the coordinates as values.
    
    Parameters:
    lenGrid (int): The number of grid points in one dimension.
    
    Returns:
    posLabels (dict): A dictionary with the labels of the grid points as keys and the coordinates as values.
    '''
    posLabels = {c: [k1, k0] for c, (k0, k1) in enumerate(((k0, k1) for k0 in range(1, lenGrid + 1) for k1 in range(1, lenGrid + 1)), 1)}
    return posLabels

def ConeMaker(turbineI, dist, r, d, lenGrid):
    '''
    Creates a cone of influence for a given turbine based on distance, radius, and direction.
    
    Parameters:
    turbineI (int): The index of the turbine.
    dist (float): The distance parameter.
    r (float): The radius parameter.
    d (list): A list containing direction and other parameters.
    lenGrid (int): The number of grid points in one dimension.
    
    Returns:
    cone (dict): A dictionary representing the cone of influence.
    '''
    theta = math.radians(d[0])
    cosTheta, sinTheta = math.cos(theta), math.sin(theta)
    pos = Labelling(lenGrid)
    ith = pos[turbineI]
    cone = {k: 0 for k in range(1, lenGrid ** 2 + 1)}
    cone[turbineI] = 1

    for k, (x, y) in pos.items():
        xNew = x * cosTheta - y * sinTheta
        yNew = x * sinTheta + y * cosTheta
        pos[k] = [xNew, yNew]
        if k != turbineI and (yNew - ith[1]) <= (dist + 0.001) and abs(xNew - ith[0]) < abs(r * (yNew - ith[1])):
            cone[k] = 1

    for k, (x, y) in pos.items():
        if y - ith[1] < 0:
            cone[k] = 0

    return cone

def Jansen(i, j, D, wfPara, lenGrid):
    '''
    Calculates the reduced wind speed at turbine j due to the wake effect from turbine i.
    
    Parameters:
    i (int): The index of the first turbine.
    j (int): The index of the second turbine.
    D (list): A list containing direction and other parameters.
    wfPara (dict): A dictionary containing wind farm parameters.
    lenGrid (int): The number of grid points in one dimension.
    
    Returns:
    reducedSpeed (float): The reduced wind speed at turbine j.
    '''
    turbineRadius = 0.33  # turbine radius is 1/3 of the site size
    pos = Labelling(lenGrid)
    I, J = pos[i], pos[j]
    dist = math.sqrt((I[0] - J[0]) ** 2 + (I[1] - J[1]) ** 2)
    freeSpeed = D[1]
    alpha = (wfPara['r'] - turbineRadius) / wfPara['x']
    a = 0.1
    reducedSpeed = freeSpeed * (1 - 2 * a / (1 + alpha * (dist / wfPara['r']) ** 2) ** 2)
    return reducedSpeed

def ReducedJansenFactor(i, j, d, wfPara, lenGrid):
    '''
    Calculates the reduced wind speed factor at turbine j due to the wake effect from turbine i.
    
    Parameters:
    i (int): The index of the first turbine.
    j (int): The index of the second turbine.
    d (list): A list containing direction and other parameters.
    wfPara (dict): A dictionary containing wind farm parameters.
    lenGrid (int): The number of grid points in one dimension.
    
    Returns:
    RF (float): The reduced wind speed factor at turbine j.
    '''
    cone = ConeMaker(i, wfPara['x'], wfPara['r'], d, lenGrid)
    if cone[j] == 1:
        return Jansen(i, j, d, wfPara, lenGrid)
    return d[1]

def ConeMatrix(cone, lenGrid):
    '''
    Converts the cone of influence into a matrix representation.
    
    Parameters:
    cone (dict): A dictionary representing the cone of influence.
    lenGrid (int): The number of grid points in one dimension.
    
    Returns:
    M (numpy.ndarray): A matrix representation of the cone of influence.
    '''
    M = np.zeros((lenGrid, lenGrid))
    for k, v in cone.items():
        k0, k1 = divmod(k - 1, lenGrid)
        M[k0, k1] = v
    return M

def NumberConstraint(m, lenGrid):
    '''
    Creates a number constraint matrix for the QUBO problem.
    
    Parameters:
    m (int): The maximum number of turbines.
    lenGrid (int): The number of grid points in one dimension.
    
    Returns:
    C (numpy.ndarray): The number constraint matrix.
    '''
    n = lenGrid ** 2
    C = np.zeros((n, n))
    
    for i in range(n):
        C[i, i] = 1
    for i in range(n):
        for j in range(i+1, n):
            C[i, j] = C[j, i] = 1 

    penalty = 2 * m
    
    return C - penalty * np.eye(n)

def ProximityConstraint(E, lenGrid):
    '''
    Creates a proximity constraint matrix for the QUBO problem.
    
    Parameters:
    E (float): The proximity threshold.
    lenGrid (int): The number of grid points in one dimension.
    
    Returns:
    C (numpy.ndarray): The proximity constraint matrix.
    '''
    C = np.zeros((lenGrid ** 2, lenGrid ** 2))
    labels = Labelling(lenGrid)
    for i in range(lenGrid ** 2):
        for j in range(lenGrid ** 2):
            if i != j:
                x1, x2 = labels[i + 1], labels[j + 1]
                distance = math.sqrt((x1[0] - x2[0]) ** 2 + (x1[1] - x2[1]) ** 2)
                if distance <= E:
                    C[i, j] = 1
    return C

def WindfarmQ(parameters):
    '''
    Creates the QUBO matrix for the wind farm layout optimization problem.
    
    Parameters:
    parameters (dict): A dictionary containing wind farm parameters.
    
    Returns:
    Q (numpy.ndarray): The QUBO matrix.
    '''
    D = parameters['D']
    x, r, m, E = parameters['x'], parameters['r'], parameters['m'], parameters['E']
    lam1, lam2 = parameters['lam1'], parameters['lam2']
    lenGrid = parameters['len_grid']
    Q = np.zeros((lenGrid ** 2, lenGrid ** 2))
    
    for d in D:  
        for k in range(1, lenGrid ** 2 + 1):
            for j in range(1, lenGrid ** 2 + 1):
                if k != j:
                    reducedSpeed = ReducedJansenFactor(k, j, d, parameters, lenGrid)
                    term = d[2] ** 3 - reducedSpeed ** 3
                    Q[k - 1, j - 1] += -1 / 3 * d[2] * term
                else:
                    Q[k - 1, j - 1] += 1 / 3 * d[2] * (d[1] ** 3)
    
    f = 0.5 * (Q + Q.T)
    g = lam1*NumberConstraint(m, lenGrid) + lam2*ProximityConstraint(E, lenGrid)
    C = -f + g
    return C

def Energy(parameters, solution):
    '''
    Calculates the energy of a given wind farm layout solution.
    
    Parameters:
    parameters (dict): A dictionary containing wind farm parameters.
    solution (list): A list representing the wind farm layout solution.
    
    Returns:
    energy (float): The energy of the given solution.
    '''
    D = parameters['D']
    x, r, E = parameters['x'], parameters['r'], parameters['E']
    lenGrid = parameters['len_grid']
    energy = 0
    
    for d in D:
        freeSpeed, prob = d[1], d[2]
        for k in range(1, len(solution) + 1):
            if solution[k - 1] == 1:
                cone = ConeMaker(k, x, r, d, lenGrid)
                term = sum(freeSpeed ** 3 - ReducedJansenFactor(k, c, d, parameters, lenGrid) ** 3 for c in cone if cone[c] == 1 and solution[c - 1] and c != k)
                energy += 1 / 3 * prob * (freeSpeed ** 3 - term)
    return energy

def WakeSpeeds(parameters, turbineI, d):
    '''
    Calculates the wake speeds for a given turbine.
    
    Parameters:
    parameters (dict): A dictionary containing wind farm parameters.
    turbineI (int): The index of the turbine.
    d (list): A list containing direction and other parameters.
    
    Returns:
    mat2 (numpy.ndarray): A matrix representing the wake speeds.
    '''
    dist, r = parameters['x'], parameters['r']
    lenGrid = parameters['len_grid']
    cone = ConeMaker(turbineI, dist, r, d, lenGrid)
    mat = ConeMatrix(cone, lenGrid)
    mat2 = np.full((lenGrid, lenGrid), d[1])
    
    for u in range(lenGrid):
        for v in range(lenGrid):
            if mat[v, u] == 1:
                c = u * lenGrid + v + 1
                mat2[v, u] = d[1] if c == turbineI else ReducedJansenFactor(turbineI, c, d, parameters, lenGrid)
    
    return mat2
