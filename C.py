#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jun 19 12:02:47 2021

@author: shaktimukker
"""

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

rhoB=epi*rhoA+(1-epi)*rhoS
Kb=epi*Ka+(1-epi)*Ks

CCP0=0