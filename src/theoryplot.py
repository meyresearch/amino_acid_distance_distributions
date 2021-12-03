#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 26 10:18:53 2021
@author: nora
"""
import networkx as nx
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import random
from scipy.optimize import curve_fit

def nthHarmonic(n) :
    N=int(n)
    harmonic = 1.00
    for i in range(2, N + 1) :
        harmonic += 1 / i
 
    return harmonic

def fhigh(x, N, HN2):
     #return 2/N/(HN2-1)*(x-1) - 2/N #setting the sum to its max
     #return (2/N/(HN2-1)+2/N)*(x-1) #setting the sum to its min
     #return (2/N/(HN2-1)+1/N)*(x-1) -1/N #setting the sum to mean
     return 2/N *((HN2 - nthHarmonic(x) + 1)/(HN2-1)*x - HN2/(HN2-1)) # using the actual sum
 
def flow(x, N):
     return -2/N**2*x**2+(2/N+6/N**2)*x-(2/N+4/N**2)
 
def P(x, N, HN2, a, A):
     return ((1-flow(x, N))/(1-fhigh(x, N, HN2)))**a *A/fhigh(x, N, HN2)

def Ptest(x, N, HN2, a, A):
     return A/fhigh(x, N, HN2)
 
def Powerlaw(x):
     return N*1/pow(x,1.2)
 
def Pin(x,N,HN2,a,A,k):
    C=sum([((1-flow(x, N))/(1-fhigh(x, N, HN2)))**a*(1-fhigh(x, N, HN2))**k for x in sumrange])
    return ((1-flow(x, N))/(1-fhigh(x, N, HN2)))**a*(1-fhigh(x, N, HN2))**k/C

def antiprop(x,con,alpha):
    return con/x/(HN2-1)

#setting all the constants. N is the chain length, HN2 is the harmonic number of N/2, because it appears many times. A is the factor from the geometric series, a is how many steps we use Pin from eq.(7)
#rather than antiprop, sumrange is the range we alrays sum over, plotrange is a shorter range we use for the last 
#plot so the axes don't get scaled such that we don't see anything any more
N=1000
HN2=nthHarmonic(N//2)
A=1
a=6
sumrange=range(2,N//2)
plotrange=range(20,N//2)

#plot of flow and fhigh respectively to check that they are between 0 and 1, they used to not be, because of all my mistakes
plt.plot(sumrange,[fhigh(s,N,HN2) for s in sumrange],label='f high')
plt.plot(range(1,N//2),flow(np.array(range(1,N//2)),N),label='f low')
plt.legend(loc=1,prop={'size':12})

#This is the important figure of the theoretic result. Specifically the blue line is. the orange line is the result for a=0,
#so a more simplified approximation. The green line is a power law with exponent 1.2 because i vaguely remembered that was 
#what we found in the data, so i decided to use that for comparison.
plt.figure()
plt.yscale('log')
plt.xscale('log')
plt.plot(sumrange,[P(s,N,HN2,a,A) for s in sumrange],label='P')
#plt.plot(plotrange,P(np.array(plotrange),N,HN2,a,A),label='P')
plt.plot(sumrange,[Ptest(s,N,HN2,a,A) for s in sumrange],label='Ptest')
plt.plot(sumrange,Powerlaw(np.array(sumrange)),label='powerlaw')
plt.legend(loc=1,prop={'size':12})


#Here i plotted the approximation of Pin from equation (8) against the resulting Pin's for different values of k to see when 
#and how it converges. It appears that it converges not that convincingly/only on average (that one i haven't checked but since
#a mean involves a sum over all k and we end up with a sort of power law after the sum in eq.12 that seems reasonable). I think
#I would just take this as a not great but look the end result is fine sort of approximation. We can put a figure like this in 
#the supplement if we want to or wait until a reviewer makes us. In terms of computational complexity i find it hard to imagine
#that an approximation any more complex would still be viable for analytical calculations.
plt.figure()
#plt.yscale('log')
#plt.xscale('log')
plt.plot(plotrange,[2/N for s in plotrange],label='k=0')
plt.plot(plotrange,[Pin(s,N,HN2,a,A,2) for s in plotrange],label='k=2')
plt.plot(plotrange,[Pin(s,N,HN2,a,A,4) for s in plotrange],label='k=4')
plt.plot(plotrange,[Pin(s,N,HN2,a,A,8) for s in plotrange],label='k=8')
plt.plot(plotrange,[Pin(s,N,HN2,a,A,32) for s in plotrange],label='k=32')
plt.plot(plotrange,[antiprop(s,1,1) for s in plotrange],label='1/s')

plt.legend(loc=1,prop={'size':12})


plt.show()