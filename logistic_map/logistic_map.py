#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 26 19:26:04 2022

@author: williamhutzel
"""

import matplotlib.pyplot as plt
import numpy as np
import itertools

def logistic_map(mu,n_start,iterations):
    population = [n_start]
    time = [0]
    for i in range(1,iterations):
        population.append(mu*population[i-1]*(1-population[i-1]))
        time.append(i)
    return time,population

def logistic_map_last_gens(mu,n_start,iterations):
    population = [n_start]
    generations = []
    for i in range(1,iterations):
        population.append(round(mu*population[i-1]*(1-population[i-1]),4))
        if (i> 200) and (population[i]>0):
            generations.append([round(mu,4),population[i]])
    generations = list(generations for generations,_ in itertools.groupby(generations))
    return generations

def diff_plot():
    eps = 2e-14
    time,population = logistic_map(4, .75 , 500)
    time2,population2 = logistic_map(4*(1+eps), .75, 500)
    diff = np.subtract(population2, population)
    plt.plot(time,diff)
    
if __name__ == "__main__":
    mu_steps = 1000
    x0_steps = 100
    generations = 500
    mu_lower = 1
    mu_upper = 4
    x0_lower = .001
    x0_upper = 1
    populations = []
    
    for mu in np.linspace(mu_lower,mu_upper,mu_steps):
        for x0 in np.linspace(x0_lower,x0_upper,x0_steps):
            populations.append(logistic_map_last_gens(mu, x0 , generations))
    
    temp = [ele[0] for ele in populations if ele != []]
    temp = list(temp for temp,_ in itertools.groupby(temp))
    x = [item[0] for item in temp]
    y = [item[1] for item in temp]        
    plt.scatter(x,y,s=.01)
            
        