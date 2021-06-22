#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jun 19 12:02:47 2021

@author: shaktimukker
"""
"""
a,b,c,d=18.30, 3816.44, 227.02, 133.322
Cpa,Cps=1180,2500
IDS=.34
Ka,Ks=74.16, 1440
p,R=101325, .03
T0,Tj,Ta,Vz,W0=27, 30, 30, 60, 2.33*1000
X0,Y0,YQ,Z=2.1, 700, 8.4*1000, .6
epi,lmbda,rhoA,rhoS=.38, 2414.3, 1.14, 750 
Vz,V=60,.003
h=200

rhoB=epi*rhoA+(1-epi)*rhoS
Kb=epi*Ka+(1-epi)*Ks

CCP0=0

# parameters to be used
Xm,Um=183,.33
Yxco2,Mco2=3.3,.01
Yxw,Mw=3.3,.01
"""

a,b,c,d=18.30, 3816.44, 227.02, 133.322
Cpa,Cps=1180,4349.8
IDS=5
Ka,Ks=74.16, 219.6
p,R=101325, .2667
T0,Tj,Ta,Vz,W0=33.5, 33.5, 33.5, 6114.65, 2.33*1000
X0,Y0,YQ,Z=10, 700, 8.4*1000, 1
epi,lmbda,rhoA,rhoS=.4016, 2414.3, 1.14, 1250 
Vz,V=6114.65,19.63*(10**-8)
h=16000

rhoB=epi*rhoA+(1-epi)*rhoS
Kb=epi*Ka+(1-epi)*Ks

CCP0=0

# parameters to be used
Xm,Um=18.7,.1986
Yxco2,Mco2=5.585,.00554989
Yxw,Mw=3.3,.01

