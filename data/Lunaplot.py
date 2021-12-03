import networkx as nx
import numpy as np
import scipy as sp
import scipy.stats
import matplotlib.pyplot as plt
import matplotlib
import math as m
import random
import time
import bisect
from scipy.optimize import curve_fit
#from compiler.ast import flatten
from operator import add

def func(x, a, b):
     return a * x**(-b)

def func1(x, a):
     return a/x

def expfit(x, a, b):
     return a * np.exp(x*(-b))

n=1000
Luna=np.zeros((10,498))
Luna[1]=np.loadtxt('luna_noencircle1')
Luna[2]=np.loadtxt('luna_noencircle2')
Luna[3]=np.loadtxt('luna_noencircle3')
Luna[4]=np.loadtxt('luna_noencircle4')
Luna[5]=np.loadtxt('luna_noencircle5')
Luna[6]=np.loadtxt('luna_noencircle6')
Luna[7]=np.loadtxt('luna_noencircle7')
Luna[8]=np.loadtxt('luna_noencircle8')
Luna[9]=np.loadtxt('luna_noencircle9')
Luna[0]=np.loadtxt('luna_noencircle10')

#n=2000
#Luna[1]=np.loadtxt('luna1')
#Luna[2]=np.loadtxt('luna2')
#Luna[3]=np.loadtxt('luna3')
#Luna[4]=np.loadtxt('luna4')
#Luna[5]=np.loadtxt('luna5')
#Luna[6]=np.loadtxt('luna6')
#Luna[7]=np.loadtxt('luna7')
#Luna[8]=np.loadtxt('luna8')
#Luna[9]=np.loadtxt('luna9')
#Luna[10]=np.loadtxt('luna10')
#Luna[11]=np.loadtxt('luna11')
#Luna[12]=np.loadtxt('luna12')
#Luna[13]=np.loadtxt('luna13')
#Luna[14]=np.loadtxt('luna14')
#Luna[15]=np.loadtxt('luna15')
#Luna[16]=np.loadtxt('luna16')
#Luna[17]=np.loadtxt('luna17')
#Luna[18]=np.loadtxt('luna18')
#Luna[19]=np.loadtxt('luna19')
#Luna[0]=np.loadtxt('luna10')

Lunatot=np.mean(Luna,axis=0)
#print(Luna[1])

plotx=np.logspace(0.3, 2.7, num=20)
#plotx=np.linspace(1, 1000, num=10)
bin_means, bin_edges, binnumber = sp.stats.binned_statistic(np.array(range(n//2-2)),np.array(Lunatot), statistic=np.mean, bins=plotx)
bin_centers = bin_edges[1:] - (bin_edges[1:] - bin_edges[0:19])/2
print(bin_means)
print(bin_centers)
#fit1, fitcov1=curve_fit(func, range(2,n/2),Luna)
#fit1, fitcov1=curve_fit(func, bin_edges[1:20],bin_means[0:20])
fit1, fitcov1=curve_fit(func, bin_centers[3:], bin_means[3:])

print(fit1)

#Plot linklength distribution
#plt.plot(Luna[0],color='lightgrey', marker='.',linewidth=0)
#plt.plot(Luna[1],color='lightgrey', marker='.',linewidth=0)
#plt.plot(Luna[2],color='lightgrey', marker='.',linewidth=0)
#plt.plot(Luna[3],color='lightgrey', marker='.',linewidth=0)
#plt.plot(Luna[4],color='lightgrey', marker='.',linewidth=0)
#plt.plot(Luna[5],color='lightgrey', marker='.',linewidth=0)
#plt.plot(Luna[6],color='lightgrey', marker='.',linewidth=0)
#plt.plot(Luna[7],color='lightgrey', marker='.',linewidth=0)
#plt.plot(Luna[8],color='lightgrey', marker='.',linewidth=0)
#plt.plot(Luna[9],color='lightgrey', marker='.',linewidth=0)
#plt.plot(Luna[10],color='lightgrey', marker='.',linewidth=0)
#plt.plot(Luna[11],color='lightgrey', marker='.',linewidth=0)
#plt.plot(Luna[12],color='lightgrey', marker='.',linewidth=0)
#plt.plot(Luna[13],color='lightgrey', marker='.',linewidth=0)
#plt.plot(Luna[14],color='lightgrey', marker='.',linewidth=0)
#plt.plot(Luna[15],color='lightgrey', marker='.',linewidth=0)
#plt.plot(Luna[16],color='lightgrey', marker='.',linewidth=0)
#plt.plot(Luna[17],color='lightgrey', marker='.',linewidth=0)
#plt.plot(Luna[18],color='lightgrey', marker='.',linewidth=0)
#plt.plot(Luna[19],color='lightgrey', marker='.',linewidth=0)

##Plot linklength distribution
plt.plot(Luna[0],color='lightgrey')
plt.plot(Luna[1],color='lightgrey')
plt.plot(Luna[2],color='lightgrey')
plt.plot(Luna[3],color='lightgrey')
plt.plot(Luna[4],color='lightgrey')
plt.plot(Luna[5],color='lightgrey')
plt.plot(Luna[6],color='lightgrey')
plt.plot(Luna[7],color='lightgrey')
plt.plot(Luna[8],color='lightgrey')
plt.plot(Luna[9],color='lightgrey')
#plt.plot(Luna[10],color='lightgrey')
#plt.plot(Luna[11],color='lightgrey')
#plt.plot(Luna[12],color='lightgrey')
#plt.plot(Luna[13],color='lightgrey')
#plt.plot(Luna[14],color='lightgrey')
#plt.plot(Luna[15],color='lightgrey')
#plt.plot(Luna[16],color='lightgrey')
#plt.plot(Luna[17],color='lightgrey')
#plt.plot(Luna[18],color='lightgrey')
#plt.plot(Luna[19],color='lightgrey')

plt.plot(Lunatot)

plt.plot(range(1,n//2),func1(np.array(range(1,n//2)),Lunatot[1]),label='analytics...',color='orange',linewidth=2, linestyle='--')
plt.plot(range(1,n//2),func(np.array(range(1,n//2)),fit1[0],fit1[1]),label='fit: a=%.2f, b=%.2f'%(fit1[0],fit1[1]),color='red',linewidth=2)
#plt.errorbar(bin_edges[0:9],bin_means, marker='o',linewidth=0,label='data',color='blue',elinewidth=1,barsabove=True)
#plt.errorbar(bin_edges[1:10],bin_means, marker='o',linewidth=0,label='data',color='green',elinewidth=1,barsabove=True)
plt.errorbar(bin_centers,bin_means, marker='o',linewidth=0,label='data',color='black',elinewidth=1,barsabove=True)
plt.xlabel('Link length', fontsize=20)
plt.ylabel('Frequency', fontsize=20)
plt.tight_layout()
plt.xlim([1,1000])
plt.ylim([0.1,500])
plt.yscale('log')
plt.xscale('log')
plt.legend(loc=3,prop={'size':12})
#plt.savefig('lilinew_lines.pdf')

plt.show()