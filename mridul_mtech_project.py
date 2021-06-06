# -*- coding: utf-8 -*-
"""
Created on Sat May 29 16:31:46 2021

@author: shakt
"""
import odenlls as nlr
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import scipy as sp

path=r"C:\Users\shakt\OneDrive\Documents\MyProjects\projects\mridul_mtech_project\\"
data=pd.read_excel(path+"data.xlsx")

#1
# fit a second degree polynomial to the economic data
from numpy import arange
from pandas import read_csv
from scipy.optimize import curve_fit
from matplotlib import pyplot
from gekko import GEKKO
    
# eq8
def fn_eq8(X,a1,a2,b0,b1,b2):
    Y = X[:, 0]
    T = X[:, 1]
    Tmin,Tmax=18,45
    #a1,a2,b0,b1,b2=para_eq8[0],para_eq8[1],para_eq8[2],para_eq8[3],para_eq8[4]
    f=((a1*(T-Tmin)*(1-np.exp(a2*(T-Tmax))))**2)*(b0+b1*Y+b2*(Y**2))
    return f
#eq9
def fn_eq9(X,a1,a2,b0,b1,b2):
    Y = X[:, 0]
    T = X[:, 1]
    Tmin,Tmax=18,45
    #a1,a2,b0,b1,b2=para_eq8[0],para_eq8[1],para_eq8[2],para_eq8[3],para_eq8[4]
    f=(a1*(T-Tmin)*(1-np.exp(a2*(T-Tmax)))**2)*(b0+b1*Y+b2*(Y**2))
    return f
    

guesses = [2.981,0.201,0.431,2.010,0.528]
guesses=[2.189, 0.284, 0.418, 1.220 ,0.867]
X=np.array(data[["Y","T"]].values.tolist())
f=data["Um"].values.tolist()
popt, pcov = curve_fit(fn_eq9, X, f, guesses)
sigma = np.sqrt(np.diag(pcov))

#plot 
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt

X = np.column_stack(X[:,0], X[:,1]) # independent variables

f = fn_eq8(X,*popt)

fig = plt.figure()
ax = fig.gca(projection = '3d')

ax.plot(X[:,0], X[:,1], f)
ax.set_xlabel('Y')
ax.set_ylabel('T')
ax.set_zlabel('f(Y,T)')

plt.savefig('images/graphical-mulvar-1.png')



