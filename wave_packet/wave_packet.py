#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 10 19:58:47 2022

@author: williamhutzel
"""
import math
import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import UnivariateSpline
from scipy.fft import fft, fftfreq


def y1(t):
    return math.exp(-(t*t)/2)

def y2(t):
    return math.sin(8*t)*math.exp(-(t*t)/2)

def y3(t):
    return (1-t*t)*math.exp(-(t*t)/2)

def fwhm(time,f_values):
    spline = UnivariateSpline(time, f_values-.5*np.max(f_values),s=0)
    r1,r2 = spline.roots()
    return r1,r2

def runy1():
    """
    This needs to be refactored to compensate for the bit flipping
    done using the fast fourier transform
    Returns
    -------
    None.

    """
    y1_val = []
    t = np.linspace(-10,10,num=100)
    freqs_arranged = np.zeros(100)
    vals_arranged =  np.zeros(100)
    
    
    for time in t:
        y1_val.append(y1(time))
        
    plt.figure(0)
    plt.plot(t,y1_val)
    r1,r2 = fwhm(t,y1_val)
    print('r1(t): ',r1,'\t r2(t): ',r2)
    sp = np.fft.fft(y1_val)
    
    plt.figure(1)
    xf = fftfreq(100,1/5)
    freq = np.fft.fftfreq(t.shape[-1])
    yf = fft(y1_val)
    plt.plot(xf,sp.real,freq,sp.imag)
    
    vals_arranged[:50] = np.abs(yf[50:])
    vals_arranged[50:] = np.abs(yf[:50])
    freqs_arranged[0:50] = xf[50:]
    freqs_arranged[50:] = xf[:50]
    
    plt.figure(2)
    #abs_y1_val = np.abs(sp)
    plt.plot(freqs_arranged,vals_arranged)
    r1,r2 = fwhm(freqs_arranged,vals_arranged)
    print('r1(w): ',r1,'\t r2(w): ',r2)
    
def runy2():
    y1_val = []
    t = np.linspace(-10,10,num=100)
    
    
    for time in t:
        y1_val.append(y2(time))
        
    plt.plot(t,y1_val)
    
def runy3():
    y1_val = []
    t = np.linspace(-10,10,num=100)
    
    
    for time in t:
        y1_val.append(y3(time))
        
    plt.plot(t,y1_val)


if __name__ == "__main__":
    runy1()
    # y1_val = []
    # t = np.linspace(-10,10,num=100)
    
    
    # for time in t:
    #     y1_val.append(y1(time))
        
    # plt.figure(0)
    # plt.plot(t,y1_val)
    # r1,r2 = fwhm(t,y1_val)
    # print('r1(t): ',r1,'\t r2(t): ',r2)
    # sp = np.fft.fft(y1_val)
    
    # plt.figure(1)
    # freq = np.fft.fftfreq(t.shape[-1])
    # plt.plot(freq,sp.real,freq,sp.imag)
    
    # plt.figure(2)
    # abs_y1_val = np.sqrt(sp.real*sp.real+sp. imag*sp.imag)
    # plt.plot(freq,abs_y1_val)
    # r1,r2 = fwhm(freq,abs_y1_val)
    # print('r1(w): ',r1,'\t r2(w): ',r2)