import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

def resl_matrix(ts,N=10):
    '''input: ts and N, output: var,cij'''
    # 线性常微分方程组的矩阵解法
    sigma=0.5**0.5
    x0=np.array([[1],[1]])
    matrix=np.array([[-2*(1-sigma**2),2],[2/(N-1),-2/(N-1)]])
    lam,um=np.linalg.eig(matrix)
    uminv=np.linalg.inv(um)
    t=0

    gam=np.array([[np.exp(lam[0]*t),0],[0,np.exp(lam[1]*t)]])
    resol=um@(gam)@uminv@x0
    y=[]
    for t in ts:
        gam=np.array([[np.exp(lam[0]*t),0],[0,np.exp(lam[1]*t)]],)
        resol=um@(gam)@uminv@x0
        y.append(resol.flatten())
    ym=np.asarray(y)
    return ym[:,0]-1,np.divide(ym[:,1]-1,ym[:,0]-1)



def gxyt(y,t,N=10):
    sigma_2=0.5
    f1 = 2*(y[1]-(1-sigma_2)*y[0])
    f2 = 2*(y[0]-y[1])/(N-1)
    return np.array([f1,f2])

def rel_rk(tau,T=10000,N=10):
    '''input: dt, T and N, output: ts, var,cij'''
    # bushu = peri_t*20/tau
    bushu = T/tau
    t = [0,]
    y=np.array([1,1])
    # yn = [y[0],]
    yn = [y,]
    for i in range(int(bushu)):
        c11 = gxyt(y,t[-1],N=N)*tau
        c12 = gxyt(y + c11*0.5,t[-1] + tau*0.5,N=N)*tau
        c13 = gxyt(y + c12*0.5,t[-1] + tau*0.5,N=N)*tau
        c14 = gxyt(y + c13,t[-1] + tau,N=N)*tau
        dy1 = (c11+2*c12+2*c13+c14)/6.

        y = y + dy1
        t.append(t[-1] + tau)
        # yn.append(y[0])
        yn.append(y)
    ym = np.asarray(yn)
    return t,ym[:,0]-1, np.divide((ym[:,1]-1),ym[:,0]-1)

def gxyt_j(y,t,N=10):
    sigma_2=0.5
    f1 = 2*(-y[1]-(-1-sigma_2)*y[0])
    f2 = - 2*(y[0]-y[1])/(N-1)
    return np.array([f1,f2])


def rel_rkj(tau,T=10000,N=10):
    '''input: dt, T and N, output: ts, var,cij'''
    # bushu = peri_t*20/tau
    bushu = T/tau
    t = [0,]
    y=np.array([1,1])
    # yn = [y[0],]
    yn = [y,]
    for i in range(int(bushu)):
        c11 = gxyt_j(y,t[-1],N=N)*tau
        c12 = gxyt_j(y + c11*0.5,t[-1] + tau*0.5,N=N)*tau
        c13 = gxyt_j(y + c12*0.5,t[-1] + tau*0.5,N=N)*tau
        c14 = gxyt_j(y + c13,t[-1] + tau,N=N)*tau
        dy1 = (c11+2*c12+2*c13+c14)/6.

        y = y + dy1
        t.append(t[-1] + tau)
        # yn.append(y[0])
        yn.append(y)
    ym = np.asarray(yn)
    return t,ym[:,0]-1, np.divide((ym[:,1]-1),ym[:,0]-1)
