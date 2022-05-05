#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 26 19:26:04 2022

@author: williamhutzel
"""

import matplotlib.pyplot as plt
import numpy as np

def logistic_map(mu,n_start,iterations):
    population = [n_start]
    time = [0]
    for i in range(1,iterations):
        population.append(mu*population[i-1]*(1-population[i-1]))
        time.append(i)
    return time,population
    
if __name__ == "__main__":
    eps = 2e-14
    time,population = logistic_map(4, .75 , 500)
    time2,population2 = logistic_map(4*(1+eps), .75, 500)
    diff = np.subtract(population2, population)
    plt.plot(time,diff)