"""Analytic derivatives of a two-layer tanh regressor; no autodiff dependency."""
import numpy as np


class MLP:
    def __init__(self, width):
        self.width = width
        self.npar = 4*width+1

    def initial(self, seed):
        rng = np.random.default_rng(seed)
        W = rng.normal(size=(self.width,2))
        b = rng.normal(scale=.3,size=self.width)
        a = rng.normal(scale=.5,size=self.width)
        return np.r_[np.column_stack([W,b]).ravel(),a,0.]

    def data(self, theta, X):
        h = self.width; weights = theta[:3*h].reshape(h,3); a = theta[3*h:4*h]
        xa = np.column_stack([X,np.ones(len(X))])
        t = np.tanh(xa@weights.T); s = 1-t*t
        q = t@a+theta[-1]
        J = np.column_stack([(s[:,:,None]*a[None,:,None]*xa[:,None,:]).reshape(len(X),3*h),t,np.ones(len(X))])
        return q,J,(xa,t,s,a)

    def weighted_hessian(self, residual, cache):
        xa,t,s,a = cache; h = self.width
        H = np.zeros((self.npar,self.npar),dtype=np.result_type(residual,a,t))
        for j in range(h):
            block = slice(3*j,3*j+3); aj = 3*h+j
            weights = residual*a[j]*(-2*t[:,j]*s[:,j])
            H[block,block] = xa.T@(weights[:,None]*xa)
            cross = xa.T@(residual*s[:,j])
            H[block,aj] = cross; H[aj,block] = cross
        return H


class FixedFeatures:
    def __init__(self, base, theta):
        h = base.width; self.weights = theta[:3*h].reshape(h,3).copy()
        self.npar = h+1; self.initial_theta = np.r_[theta[3*h:4*h],theta[-1]]

    def data(self, theta, X):
        xa = np.column_stack([X,np.ones(len(X))])
        features = np.column_stack([np.tanh(xa@self.weights.T),np.ones(len(X))])
        return features@theta,features,None

    def weighted_hessian(self, residual, cache):
        return np.zeros((self.npar,self.npar))


def field(model, theta, X, target, derivative=False):
    q,J,cache = model.data(theta,X); e = q-target; n = len(X)
    F = -J.T@e/n
    if derivative:
        return F,-(J.T@J+model.weighted_hessian(e,cache))/n
    return F


def inputs(config):
    a = np.arange(3)*2*np.pi/3
    train = config['input_train_radius']*np.column_stack([np.cos(a),np.sin(a)])
    control = config['input_control_radius']*np.column_stack([np.cos(a+np.pi/3),np.sin(a+np.pi/3)])
    directions = {'u':np.array([1.,-1.,0.])/np.sqrt(2),
                  'v':np.array([1.,1.,-2.])/np.sqrt(6),
                  'z':np.ones(3)/np.sqrt(3)}
    return train,np.vstack([train,control]),directions
