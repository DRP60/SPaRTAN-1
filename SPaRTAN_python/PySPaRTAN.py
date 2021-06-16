        
# -*- coding: utf-8 -*-
"""
Created on Sat Feb 29 11:44:06 2020

@author: XIM33
"""
import numpy as np
import cythKronPlus as krnP
import cythLeastR as leastR
import scipy.linalg
import functools
import time
import gc
import matplotlib.pyplot as plt

def catch_exception(f):
    @functools.wraps(f)
    def func(*args, **kwargs):
        try:
            return f(*args, **kwargs)
        except Exception as e:
            print ('Caught an exception in', f.__name__)
    return func

class ErrorCatcher(type):
    def __new__(cls, name, bases, dct):
        for m in dct:
            if hasattr(dct[m], '__call__'):
                dct[m] = catch_exception(dct[m])
        return type.__new__(cls, name, bases, dct)


class Timer:
    def __init__(self):
        self.start = time.time()

    def restart(self):
        self.start = time.time()

    def get_time_hhmmss(self):
        end = time.time()
        m, s = divmod(end - self.start, 60)
        h, m = divmod(m, 60)
        time_str = "%02d:%02d:%02d" % (h, m, s)
        return time_str

def normalize_column(A, T=0):
    if (T== 0):
        return np.divide(A,np.sqrt(np.sum(A**2, 0)))
    else:
        At = np.transpose(A)
        return np.transpose(np.divide(At,np.sqrt(np.sum(At**2, 0)))) 

class PySPaRTAN:
    __metaclass__ = ErrorCatcher
    # def __init__(self)

    def fit(self,D, P_train, Y_train, lamda = 0.1, rsL2 = 0, spectrumA = 0.95, spectrumB = 0.9 ):

        self.D = D
        self.P_train = P_train
        self.Y_train = Y_train
 
       # transformation
        A = self.Y_train.T @ self.D
        B = self.P_train.T
        Y = self.Y_train.T @ self.Y_train
    
        # SVD(A) SVD(B)  
        UA, SA, VhA = np.linalg.svd(A); VA = VhA.T
        UB, SB, VhB = np.linalg.svd(B); VB = VhB.T  
    
        a_cum_spectrum = np.cumsum(SA) / sum(SA)
        b_cum_spectrum = np.cumsum(SB) / sum(SB)
    
        da = np.nonzero(a_cum_spectrum >= spectrumA)[0][0] + 1  
        db = np.nonzero(b_cum_spectrum >= spectrumB)[0][0] + 1     
        
        Ua = UA[:, :da]
        Sa = SA[:da]; 
        Va = VA[:, :da]

        Ub = UB[:, :db]
        Sb = SB[:db]; 
        Vb = VB[:, :db]

        Yv = (Y.T).flatten()
        
        Vb = Vb.copy(order='C')
        Ua = Ua.copy(order='C')
        L = krnP.kron(Vb, Ua)

        d = np.eye(Y.shape[0],Y.shape[1])     
        cidex = np.where(d.flatten()!=0)
        diag = np.array(cidex,dtype=np.int32).flatten()
        
        Yv = Yv.copy(order='C') #make it c-like contiguous array
        diag = diag.copy(order='C')
        
        L,Yv = krnP.removeDiagC(L,Yv,diag)

        opts = dict()   
        opts['rsL2']=0

        #reshape Yv to 2darry
        Yv = Yv.reshape(Yv.shape[0],1)
        beta, b = leastR.LeastR(L,Yv,lamda, opts)        
        
        del L, Yv
        gc.collect()           
        
        self.beta = beta
        self.Ua = Ua
        self.Ub = Ub
        self.Sa = np.diag(Sa)
        self.Sb = np.diag(Sb)
        self.Va = Va
        self.Vb = Vb
        self.lamda = lamda

    def ar_model2w(self):
        m1 = self.Va
        m2 = np.linalg.pinv(self.Sa)
        m3 = self.beta.reshape(self.Va.shape[1], self.Ub.shape[1], order="F" )
        m4 = np.linalg.pinv(self.Sb)
        m5 = self.Ub.T
        ww = m1 @ m2 @ m3 @ m4 @ m5 
        return ww
    
    def ar_reconstruction(self, pred_test=None):
        A = self.Y_train.T @ pred_test
        O = scipy.linalg.orth(self.Y_train)
        cm = scipy.linalg.lstsq(O, self.Y_train)[0]
        ct = scipy.linalg.lstsq(cm.T, A)[0]
        pred = O @ ct
        return pred
        
        
    def predict(self, P_test=None):
        
        if P_test is not None:
            self.P_test = P_test
                       
        w = self.ar_model2w()
        pred = self.D @ (w @ self.P_test.T)
           
        aff_rec = self.ar_reconstruction(pred)
        
        self.Y_pred = aff_rec
        return aff_rec
    
    # get the correlation between predicted Y and its o
    def get_corr(self,Y_pred, Y_test, plot=False):
             
        corr = np.corrcoef(Y_test.ravel(order='F'), Y_pred.ravel(order='F'))[0,1]
        if plot:
            plt.plot(Y_test.ravel(order='F'), Y_pred.ravel(order='F'), linestyle='none', marker='+')
            plt.title('reconstruction of Y test, corr={:.2f}'.format(corr))
        
        return corr
    
    #get interaction between TF and Protein
    def get_W(self):
        self.W = self.ar_model2w()
        return self.W

    # get projected protein expression on Y, default is Y_train
    def get_projP(self, Y=None):
        if Y == None:
            Y = self.Y_train
        return (Y.T @ self.D @ self.W)	

    # get projected TF activity on P, default is P_train
    def get_projD(self, P=None):
        if P == None:
            P = self.P_train
        return self.W @ P.T	


