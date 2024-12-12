WFLO
All the code for WFLO (Jansen wake model) using quantum computing. The methods included are VQE-bases, and PCE-bases, as well as an exhaustive check. 

==============================================================================================

This project will require the following packages to run:
- numpy
- scipy
- skopt
- sklearn
- qiskit
- qiskit_aer
- itertools
- matplotlib

==============================================================================================

Project runs from main.py, controlled by a set of parameters. Below is a description of these parameters.

'len_grid': Integer
- Length of one side of the square windfarm grid.

'D': List(List(Float, Float, Float))
- Wind regime.
- Inner list format: [Angle, Speed, Probability].

'x': Float
- Length back wake will cover.
'r': Float
- Radius of wake cone per unit back.

'm': Integer
- Maximum number of turbines.

'E': FLoat   
- Proximity threshold.

'lam1': Float
- Weight for number constraint.
- Recommended: 1e4.

'lam2': Float
- Weight for proximity constraint.
- Recommended: 3e3.

'shots': Integer
- Circuit calls per quantum expected value calculation.
- Recommended: 3000.

'cvarAlpha': Float
- Proportion of samples considered when using CVaR ('VQE' only).
- Options: 0.0 to 1.0.
- Recommended: 0.8.

'solver': String
- Which quantum mapping to use.
- Options: 'PCE' or 'VQE'.

'method': String ('COBYLA' or 'Bayesian')
- Which classical optimization routine to use.
- Options: 'COBYLA' or 'Bayesian'.

'nSamples': Integer
- Number of full runs to complete

'talpha': Integer
- Alpha hyperparameter in PCE (https://arxiv.org/abs/2401.09421)
- Recommended: 10

'L': Integer
- Length of samples retained for average stopping criteria
- Recommended: 15

'tol': Float
- Tolerance for average stopping criteria
- Recommended: 1

'bayesSigma': Float
- Sigma hyperparameter in Bayesian quantum-inspired kernel (https://arxiv.org/abs/2406.06150)
- Recommended: 1000.0

'bayesGamma': Float
- Gamma hyperparameter in Bayesian quantum-inspired kernel (https://arxiv.org/abs/2406.06150)
- Recommended: 0.01

'bayesIters': Integer
- Number of Bayesian optimization steps to use.
- Recommended: 100

'ExhaustiveCheck': Boolean
- Should the exhaustive check be done.


