import os
import numpy as np

def createDir(name):
    if not os.path.isdir(name):
        os.mkdir(name)
        
def averageContiguousVals(l):
    l = np.array(l)
    l2 = np.roll(l, shift=-1)[:-1]
    l = l[:-1]
    l += (l2 - l) / 2
    return l

def EtaStr(eta):
    return str(eta).replace('.', 'p', 1)
        
