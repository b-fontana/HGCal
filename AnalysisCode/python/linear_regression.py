import json
import numpy as np
from sklearn.linear_model import LinearRegression
from CustomLinearRegressionModel import *

"""
The data is organized in 4 dimensions as follows:
1d - subdetector: 0->cee, 1->hef, 2->heb
2d - eta intevrals
3d - event
4d - event information: 0->gen energy, 1->gen eta, 2-> reco energy, 3->noise 
(see Calibration::energies_for_calibration in src/calibration.cc)
"""
def json_reader(dim1, dim2):
    data = {}
    for id1 in range(dim1):
        data[id1] = {}
        for id2 in range(dim2):
            with open('linear_regression_det'+str(id1)+'_eta'+str(id2)+'.json') as jsonfile:
                data[id1][id2] = np.array(json.load(jsonfile, parse_float=float()))
    return data

def organize_data(data):
    """puts the data in the format required by LinearRegression"""
    genen, recoen = ([] for _ in range(2))
    for id2 in range(len(data[0])):
        recoen.append( np.concatenate((data[0][id2][:,2,np.newaxis], data[1][id2][:,2,np.newaxis], data[2][id2][:,2,np.newaxis]),
                                      axis=1))
        #the genen information is repeated across subdetectors: data[1] and data[2] would also work
        genen.append(data[0][id2][:,0,np.newaxis])
    return genen, recoen

def linear_regression(x, y, custom):
    assert(len(x) == len(y))
    for i in range(len(x)):
        if custom:
            model = CustomLinearRegressionModel(X=x[i], Y=y[i])
            model.fit()
            print(model.beta)
        else:
            model = LinearRegression(fit_intercept=False, n_jobs=-1)
            model.fit(x[i], y[i])
            print(model.coef_)
        print(model.score(x[i],y[i]))

genen, recoen = organize_data(json_reader(3,1))
linear_regression(recoen, genen, False)
