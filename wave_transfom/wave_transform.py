#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr  5 20:51:16 2022

@author: williamhutzel
"""

import cmath
import numpy as np
import matplotlib.pyplot as plt

iT =  0.0;          fT =  12.0;         W = fT - iT;
N =  240;           h =  W/N
noPtsSig =  N;      noS =  20;          noTau =  90;
iTau =  0.;         iS =  0.1;          tau =  iTau;        sig =  iS 

# Need *very* small s steps for high frequency;
dTau =  W/noTau;    dS =  (W/iS)**(1./noS);
Yn =  np.zeros( (noS+1, noTau+1), float)    

def morlet(t,sig,tau):
    T = (t-tau)/sig
    val = np.exp(2j*np.pi*T*sig)*np.exp((-T*T)/2)
    return val

def mex_hat(t,sig):
   val = (1-((t*t)/(sig*(sig))))*np.exp((-t*t)/(2*sig*sig))
   return val
   
def haar_wavelet(t):
    if (t<=0 and t>-1):
        val = 1
    elif (t>0 and t<1):
        val = -1
    else:
        val = 0
    return val

def pure_sin_sig(t):
    val = np.sin(2*np.pi*t)
    return val

def sum_sin_sig(t):
    val = 2.5*np.sin(2*np.pi*t) + 6*np.sin(4*np.pi*t) + 10*np.sin(6*np.pi*t)
    return val

def nonstation_sig(t):
    if (t>=0 and t<=2):
        val = np.sin(2*np.pi*t)
    elif (t>2 and t<=8):
        val = 5*np.sin(2*np.pi*t) + 10*np.sin(4*np.pi*t)
    elif (t>8 and t<=12):
        val = 2.5*np.sin(2*np.pi*t) + 6*np.sin(4*np.pi*t) + 10*np.sin(6*np.pi*t)
    else:
        val = 0
    return val

def transform(signal,times,sig,tau):
    integral = 0
    h = max(times)/len(times)
    for i in range(0,len(signal)):
        integral += signal[i]*morlet(times[i],sig,tau)*h
    return integral / np.sqrt(sig)

def run():
    maxY =  0.001
    sig =  iS 
    time = np.linspace(-10,10,200)
    wave_val = []
    for t in time:
        wave_val.append(nonstation_sig(t))
        
    print("working, finding transform, count 20")
    for i in range( 0, noS):
        sig *= dS                                                 # Scaling
        tau = iT
        print(i)
        for j in range(0, noTau):
             tau += dTau                                      # Translate
             Yn[i, j] = transform(wave_val,time, tau, sig)
    print("transform found")  
    for i in range( 0, noS):
        for j in range( 0, noTau):
            if Yn[i, j] > maxY or Yn[i, j] < - 1 *maxY :
                maxY = abs( Yn[i, j] )                      # Find max Y       
    tau =  iT
    sig =  iS
    X = []
    Y = []
    print("normalize")      
    for i in range( 0, noS):
         Y.append(sig)
         sig *= dS                             
         for j in range( 0, noTau):
             if i == 0:
                 X.append(tau)
             tau +=   dTau                                        # Transform
             Yn[i, j] = Yn[i, j]/maxY
         tau = iT
    
        
    #plt.plot(Yn[12,:])
    x,y = np.meshgrid(X,Y)
    fig = plt.figure()
    ax = plt.axes(projection='3d')
    ax.plot_wireframe(x,y[:20],Yn[:20,:90])
    plt.contour(x,y[:90],Yn[:20,:90])      
    ax.view_init(-140, 20)
    return Yn

if __name__ == "__main__":
    trans = run()