# -*- coding: utf-8 -*-
"""
Created on Sat May 29 16:31:46 2021

@author: shakt
"""
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from gekko import GEKKO
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt

path=r"/Users/shaktimukker/Documents/MyProjects/projects/bacterial-growth-analysis/"
data=pd.read_excel(path+"data.xlsx")
    
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
X = np.column_stack(X[:,0], X[:,1]) # independent variables

f = fn_eq8(X,*popt)

fig = plt.figure()
ax = fig.gca(projection = '3d')

ax.plot(X[:,0], X[:,1], f)
ax.set_xlabel('Y')
ax.set_ylabel('T')
ax.set_zlabel('f(Y,T)')

plt.savefig('images/graphical-mulvar-1.png')


### pyomo tutorial
path=r"/Users/shaktimukker/Documents/MyProjects/projects/bacterial-growth-analysis/"
import sys
sys.path.append(path)
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.axes3d import Axes3D 
import shutil
import sys
import os.path
from pyomo.environ import *
from pyomo.dae import *
import math
import C

def model_plot(m):
    r = sorted(m.r)
    t = sorted(m.t)

    rgrid = np.zeros((len(t), len(r)))
    tgrid = np.zeros((len(t), len(r)))
    Tgrid = np.zeros((len(t), len(r)))

    for i in range(0, len(t)):
        for j in range(0, len(r)):
            rgrid[i,j] = r[j]
            tgrid[i,j] = t[i]
            Tgrid[i,j] = m.T[t[i], r[j]].value

    fig = plt.figure(figsize=(10,6))
    ax = fig.add_subplot(1, 1, 1, projection='3d')
    ax.set_xlabel('Distance r')
    ax.set_ylabel('Time t')
    ax.set_zlabel('Temperature T')
    p = ax.plot_wireframe(rgrid, tgrid, Tgrid)

m = ConcreteModel()

m.r = ContinuousSet(bounds=(0,1))
m.t = ContinuousSet(bounds=(0,2))

m.T = Var(m.t, m.r)

m.dTdt   = DerivativeVar(m.T, wrt=m.t)
m.dTdr   = DerivativeVar(m.T, wrt=m.r)
m.d2Tdr2 = DerivativeVar(m.T, wrt=(m.r, m.r))

m.pde = Constraint(m.t, m.r, rule=lambda m, t, r: m.dTdt[t,r] == m.d2Tdr2[t,r] + (1/r)*m.dTdr[t,r]
        if r > 0 and r < 1 and t > 0 else Constraint.Skip)

m.ic  = Constraint(m.r, rule=lambda m, r:    m.T[0,r] == 0)
m.bc1 = Constraint(m.t, rule=lambda m, t:    m.T[t,1] == 1 if t > 0 else Constraint.Skip)
m.bc2 = Constraint(m.t, rule=lambda m, t: m.dTdr[t,0] == 0)

TransformationFactory('dae.finite_difference').apply_to(m, nfe=20, wrt=m.r, scheme='CENTRAL')
TransformationFactory('dae.finite_difference').apply_to(m, nfe=50, wrt=m.t, scheme='BACKWARD')
SolverFactory('ipopt').solve(m).write()

model_plot(m)

##2
m = ConcreteModel()

m.r = ContinuousSet(bounds=(0,1))
m.t = ContinuousSet(bounds=(0,2))

m.T = Var(m.t, m.r)

m.dTdt   = DerivativeVar(m.T, wrt=m.t)
m.dTdr   = DerivativeVar(m.T, wrt=m.r)
m.d2Tdr2 = DerivativeVar(m.T, wrt=(m.r, m.r))

@m.Constraint(m.t, m.r)
def pde(m, t, r):
    if t == 0:
        return Constraint.Skip
    if r == 0 or r == 1:
        return Constraint.Skip
    return m.dTdt[t,r] == m.d2Tdr2[t,r]

m.obj = Objective(expr=1)

m.ic  = Constraint(m.r, rule=lambda m, r:    m.T[0,r] == 0 if r > 0 and r < 1 else Constraint.Skip)
m.bc1 = Constraint(m.t, rule=lambda m, t:    m.T[t,1] == 1)
m.bc2 = Constraint(m.t, rule=lambda m, t: m.dTdr[t,0] == 0)

TransformationFactory('dae.finite_difference').apply_to(m, nfe=50, scheme='FORWARD', wrt=m.r)
TransformationFactory('dae.finite_difference').apply_to(m, nfe=50, scheme='FORWARD', wrt=m.t)
SolverFactory('ipopt').solve(m, tee=True).write()
model_plot(m)


###curve fitting
### solve all equations together using given parameters
###solve heat equation

# parameters value for 3 liter bioreactor
a,b,c,d=18.30, 3816.44, 227.02, 133.322
Cpa,Cps=1180,2500
IDS=.34
Ka,Ks=74.16, 1440
p,R=101325, .03
T0,Tj,Ta,Vz,W0=27, 30, 30, 60, 2.33*1000
X0,Y0,YQ,Z=2.1, 700, 8.4*1000, .6
epi,lmbda,rhoA,rhoS=.38, 2414.3, 1.14, 750 

Um, Xm =.33,184
Vz=60
h=200


###solve 2D heat equation
m = ConcreteModel()

m.z = ContinuousSet(bounds=(0,C.Z))
m.r = ContinuousSet(bounds=(0,C.R))
m.t = ContinuousSet(bounds=(0,150))

m.T = Var(m.t, m.z, m.r)

m.dTdt   = DerivativeVar(m.T, wrt=m.t)
m.dTdz   = DerivativeVar(m.T, wrt=m.z)
m.d2Tdz2 = DerivativeVar(m.T, wrt=(m.z, m.z))
m.dTdr   = DerivativeVar(m.T, wrt=m.r)
m.d2Tdr2 = DerivativeVar(m.T, wrt=(m.r, m.r))


#issues
#t,T,z=0,T0,0

def X(m):
    return C.Xm/(1+((C.Xm/C.X0)-1)*math.exp(-C.Um*m.t))
def dXdt(m):
    return C.Um*X(m)*(1-X(m)/C.Xm)
def dHdT1(m):
    return .62413*C.b*C.p
def dHdT2(m):
    return (m.T+C.c)**2
def dHdT3(m):
    return C.d*np.exp(C.a-C.b/(m.T+C.c))
def dHdT4(m):
    return (C.p/dHdT3(m)-1)**2
def dHdT(m):
    return dHdT1(m)/(dHdT2(m)*dHdT4(m)*dHdT3(m))
def Cpb(m):
    return (C.epi*C.rhoA*(C.Cpa+C.lmbda*dHdT(m)+(1-C.epi)*C.rhoS*C.Cps))/C.rhoB
def heRHS(m):
    return C.rhoB*Cpb(m)*m.dTdt[t,z,r]
def heLHS1(m):
    return C.rhoS*(1-C.epi)*C.YQ*dXdt(m)
def heLHS2(m):
    return C.rhoA*C.Cpa*C.Vz*m.dTdz[t,z,r]
def heLHS3(m):
    return C.rhoA*C.lmbda*C.Vz*dHdT(m)*m.dTdz[t,z,r]
def heLHS4(m):
    return C.Kb*m.d2Tdz2[t,z,r]
def heLHS5(m):
    return C.Kb*m.dTdr[t,z,r]/m.r+C.Kb*m.d2Tdr2[t,z,r]


m.pde = Constraint(m.t, m.z, m.r, rule=lambda m, t, z, r: heRHS(m) == heLHS1(m)-heLHS2(m)-heLHS3(m)+heLHS4(m)+heLHS5(m)
        if r > 0 and r < C.R and z > 0 and z < C.Z and t > 0 else Constraint.Skip)

# =============================================================================
# m.ic  = Constraint(m.r, rule=lambda m, r:    m.T[0,r] == 0)
# m.bc1 = Constraint(m.t, rule=lambda m, t:    m.T[t,1] == 1 if t > 0 else Constraint.Skip)
# m.bc2 = Constraint(m.t, rule=lambda m, t: m.dTdr[t,0] == 0)
# 
# m.ic  = Constraint(m.z, rule=lambda m, z:    m.T[0,z] == 27)
# m.ic  = Constraint(m.r, rule=lambda m, r:    m.T[0,r] == 27)
# m.bc1 = Constraint(m.t, rule=lambda m, t:    m.dTdz[t,Z] == 0 if t > 0 else Constraint.Skip)
# m.bc2 = Constraint(m.t, rule=lambda m, t: m.T[t,0] == Ta)
# m.bc3 = Constraint(m.t, rule=lambda m, t:    m.dTdz[t,0] == 0 if t > 0 else Constraint.Skip)
# m.bc4 = Constraint(m.t, rule=lambda m, t: -Kb*m.dTdr[t,R] == h*(T-Tj))
# =============================================================================


m.ic  = Constraint(m.z,m.r, rule=lambda m, z, r:    m.T[0,z,r] == 27)
m.bc1 = Constraint(m.t,m.r, rule=lambda m, t,r:    m.dTdz[t,C.Z,r] == 0 if t > 0 else Constraint.Skip)
m.bc2 = Constraint(m.t,m.r, rule=lambda m, t,r: m.T[t,0,r] == C.Ta)
m.bc3 = Constraint(m.t,m.z, rule=lambda m, t,z:    m.dTdz[t,z,0] == 0 if t > 0 else Constraint.Skip)
m.bc4 = Constraint(m.t,m.z, rule=lambda m, t,z: -C.Kb*m.dTdr[t,z,C.R] -C.h*(m.T[t,z,C.R]-C.Tj)==0)

TransformationFactory('dae.finite_difference').apply_to(m, nfe=20, wrt=m.z, scheme='CENTRAL')
TransformationFactory('dae.finite_difference').apply_to(m, nfe=20, wrt=m.r, scheme='CENTRAL')
TransformationFactory('dae.finite_difference').apply_to(m, nfe=20, wrt=m.t, scheme='CENTRAL')
SolverFactory('ipopt').solve(m).write()

model_plot(m)

m.r
m.T[1,0,0].value

#######

###curve fitting
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from gekko import GEKKO
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt

path=r"/Users/shaktimukker/Documents/MyProjects/projects/bacterial-growth-analysis/"
data_x_vs_t=pd.read_excel(path+"data.xlsx",sheet_name="x_vs_t")
    
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
# X vs t
def x_vs_t(X,Xm,Um):
    t=X
    return Xm/(1+((Xm/C.X0)-1)*math.exp(-Um*t))

# CO2 vs t
def co2_vs_t(X,Yxco2,Mco2):
    Xm,Um=1,2
    t=X
    comn=1+((Xm/C.X0)-1)*math.exp(-Um*t)
    RHS1=1/(Yxco2*comn)
    RHS2=1/(Yxco2*(Xm/C.X0))
    RHS3=(Mco2/Um)*math.log(comn/((Xm/C.X0)*math.exp(-Um*t)))
    return C.CCP0+Xm*(RHS1-RHS2+RHS3)


    

guesses = [2.981,0.201,0.431,2.010,0.528]
guesses=[2.189, 0.284, 0.418, 1.220 ,0.867]
X=np.array(data[["Y","T"]].values.tolist())
f=data["Um"].values.tolist()
popt, pcov = curve_fit(fn_eq9, X, f, guesses)
sigma = np.sqrt(np.diag(pcov))

#plot 
X = np.column_stack(X[:,0], X[:,1]) # independent variables

f = fn_eq8(X,*popt)

fig = plt.figure()
ax = fig.gca(projection = '3d')

ax.plot(X[:,0], X[:,1], f)
ax.set_xlabel('Y')
ax.set_ylabel('T')
ax.set_zlabel('f(Y,T)')

plt.savefig('images/graphical-mulvar-1.png')

## APM
%matplotlib inline
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from scipy.optimize import curve_fit
from sklearn.metrics import r2_score

x=pd.read_csv("https://apmonitor.com/che263/uploads/Main/heart_rate.txt")
x.describe()

def bpm(t,c0,c1,c2,c3):
    return c0+c1*t-c2*np.exp(-c3*t)
g=[100,.01,100,.01]

t=x["Time (sec)"].values
hr=x["Heart Rate (BPM)"].values
g,cov=curve_fit(bpm,t,hr,g)
y =[bpm(t,g[0],g[1],g[2],g[3]) for t in x["Time (sec)"]]


plt.plot(x["Time (sec)"],x["Heart Rate (BPM)"])
plt.plot(x["Time (sec)"],y,'r.')

print('R^2: ', r2_score(y,hr))

# x vs t
%matplotlib inline
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from scipy.optimize import curve_fit
from sklearn.metrics import r2_score

path=r"/Users/shaktimukker/Documents/MyProjects/projects/bacterial-growth-analysis/"
data_x_vs_t=pd.read_excel(path+"data.xlsx",sheet_name="x_vs_t")
data_x_vs_t.describe()

# X vs t
def x_vs_t(X,Xm,Um):
    t=X
    return Xm/(1+((Xm/C.X0)-1)*np.exp(-Um*t))

g=[184,.33]
t=data_x_vs_t["t"].values
X=data_x_vs_t["x"].values
g,cov=curve_fit(x_vs_t,t,X,g)
y =[x_vs_t(t,g[0],g[1]) for t in X]

plt.plot(t,X)
plt.plot(t,y,'r')

print('R^2: ', r2_score(X,y))

####### co2 vs t
%matplotlib inline
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from scipy.optimize import curve_fit
from sklearn.metrics import r2_score
from sklearn.preprocessing import MinMaxScaler

path=r"/Users/shaktimukker/Documents/MyProjects/projects/bacterial-growth-analysis/"
data_co2_vs_t=pd.read_excel(path+"data.xlsx",sheet_name="co2_vs_t")
data_co2_vs_t.describe()

# CO2 vs t
def co2_vs_t(X,Yxco2,Mco2):
    Xm,Um=184,.33
    t=X
    comn=1+((Xm/C.X0)-1)*np.exp(-Um*t)
    RHS1=1/(Yxco2*comn)
    RHS2=1/(Yxco2*(Xm/C.X0))
    RHS3=(Mco2/Um)*np.log(comn/((Xm/C.X0)*np.exp(-Um*t)))
    return C.CCP0+Xm*(RHS1-RHS2+RHS3)

g=[2.60,.005]
t=data_co2_vs_t["t"].values
X=data_co2_vs_t["co2"].values
X=500*X
g,cov=curve_fit(co2_vs_t,t,X,g)
y =[co2_vs_t(t,g[0],g[1]) for t in X]

plt.plot(t,X)
plt.plot(t,y,'r')

print('R^2: ', r2_score(y,X))
#print('R^2: ', r2_score(X,y))







