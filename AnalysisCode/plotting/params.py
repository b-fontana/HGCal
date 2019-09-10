import numpy as np

def etavalues(samples):
    if samples == 'outer':
        return np.array([1.45, 1.49, 1.5, 1.51, 1.515, 1.52, 1.525, 1.53, 1.535, 1.54, 1.55, 1.56, 1.57, 1.58, 1.65])
    elif samples == 'inner':
        return np.array([2.7, 2.73, 2.76, 2.79, 2.82, 2.85, 2.875, 2.9, 2.91, 2.92, 2.93, 2.94, 2.975, 3.03])
    
def resvalues():
    return np.linspace(-1, 2.5, 100)
