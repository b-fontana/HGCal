import numpy as np
from scipy.optimize import minimize


def quantile_loss_generator(q):
    """generator for a q-quantile loss"""
    def aux_fun(y_true, y_pred):
        e = (y_true - y_pred)
        aux = np.where(np.greater_equal(e, 0.), q*e, (q-1)*e)
        return np.mean(aux)
    qkey='%3.2f'%q
    aux_fun.__name__ = 'quantile_loss_'+qkey.replace(".", "")
    return aux_fun


class CustomLinearRegressionModel:
    """
    Linear model with input loss function
    L2 regularization allowed
    based on https://alex.miller.im/posts/linear-model-custom-loss-function-regularization-python/
    """
    def __init__(self, X, Y, beta_init=None, 
                 loss_function=quantile_loss_generator(0.5),
                 regularization=0.0):
        
        self.loss_function = loss_function
        self.X = X
        self.Y = Y
        self.beta_init = beta_init
        self.beta = None            
        self.regularization = regularization

    def predict(self, X):
        prediction = np.matmul(X, self.beta)
        return(prediction)

    def model_error(self):
        return self.loss_function(self.predict(self.X), self.Y)
    
    def l2_regularized_loss(self, beta):
        self.beta = beta
        err=self.model_error()         
        if self.regularization==0 : return err
        l2=sum(self.regularization*np.array(self.beta)**2)
        return err+l2
            
    def fit(self, maxiter=500):        
        #initialize with dummy values
        if not self.beta_init:
            self.beta_init = np.array([1]*self.X.shape[1])
       
        if self.beta and all(self.beta_init==self.beta):
            print("Model already fit once; continuing fit with more itrations.")
            
        res = minimize(self.l2_regularized_loss, 
                       self.beta_init,
                       method='BFGS', options={'maxiter': maxiter})
        self.beta = res.x
        self.beta_init = self.beta
