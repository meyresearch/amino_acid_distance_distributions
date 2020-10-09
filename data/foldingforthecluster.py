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
from compiler.ast import flatten
from operator import add

def func(x, a, b):
     return a * x**(-b)

def expfit(x, a, b):
     return a * np.exp(x*(-b))

def pi(x):
    #return 1.0/3
    if (x<n):
        return float(4*x)/(3*(n+2*x)) +1.0/3
    else: return float(3)/2

def diameter_size(G):
    shortest_path_pairwise=nx.shortest_path(G)
    all_shortest_paths=[]

    for shortest_paths_from_src in shortest_path_pairwise.values():
        all_shortest_paths+=shortest_paths_from_src.values()

    diameter=max([len(path) for path in all_shortest_paths])
    paths_forming_diameter=[path for path in all_shortest_paths if len(path)==diameter]
    return len(set([node for path in paths_forming_diameter for node in path]))
#flattening the list and using set() to get rid of duplicate nodes that are in multiple diameters

def cluster_growing_scaling(G):
    Ml=n*np.ones((picks,int(np.sqrt(n))))
    for i in range(picks):
        n0=random.choice(range(n))
        for j in range(int(np.sqrt(n))):
            Ml[i][j]=len(nx.ego_graph(G,n0,radius=j))
    return Ml

def check(a,b,list):
    for sublist in list:
        if a in sublist and b in sublist:
            return True
    return False

def Adjust(node1,node2,A,H):
    lo=0
    for au in range(len(A)):
        if node1 in A[au] and node2 in A[au]:
            lo=A[au]
    if lo!=0:
        A.remove(lo)
        A.extend(list(nx.all_simple_paths(H.subgraph(lo),source=node1,target=node2)))
    else:
        A.extend([nx.shortest_path(H,source=node1,target=node2)])


def addlink(node1,node2,A,H,tried):
    global encircling
    if (H.degree(node1)<15 and H.degree(node2)<15) or encircling==True:
        Adjust(node1,node2,A,H)
        H.add_edge(node1,node2)
        encircling=False
    #elif H.degree(node1)==5 and H.degree(node2)<5 and encircling==False:
        #Adjust(node1,node2,A,H)
        #H.add_edge(node1,node2)
        #ring=H.neighbors(node1)
        #ringlist=[]
        #for n1 in ring:
            #for n2 in ring:
                #if check(n1,n2,A)and n1<n2:
                    #ringlist.extend([n1,n2])
                    #if H.degree(n1)<6 and H.degree(n2)<6:
                        #encircling=True
                        #addlink(n1,n2,A,H,tried)
        #outer=[x for x in ringlist if ringlist.count(x) == 1]
        #encircling=True
        #addlink(outer[0],outer[1],A,H,tried)
    #elif H.degree(node2)==5 and H.degree(node1)<5 and encircling==False:
        #Adjust(node1,node2,A,H)
        #H.add_edge(node1,node2)
        #ring=H.neighbors(node2)
        #ringlist=[]
        #for n1 in ring:
            #for n2 in ring:
                #if check(n1,n2,A)and n1<n2:
                    #ringlist.extend([n1,n2])
                    #if H.degree(n1)<6 and H.degree(n2)<6:
                        #encircling=True
                        #addlink(n1,n2,A,H,tried)
        #outer=[x for x in ringlist if ringlist.count(x) == 1]
        #encircling=True
        #addlink(outer[0],outer[1],A,H,tried)

    tried.append((node1,node2))

def Fold_my_chain(n2,gamma):
    #print  '---------------------------------------'
    starttime = time.time()
    n=int(n2)+random.randint(-3,3)
    looplens=[i*(i-3)/2 for i in range(3,n+1)]
    #create cycle graph
    H=nx.cycle_graph(n)
    tried=nx.cycle_graph(n).edges()
    alin=[]
    A=[range(n)]
    alen=[len(ii)*(len(ii)-3)/2 for ii in A]
    alin=[sum(alen[0:ii]) for ii in range(1,len(A)+1)]
    step=0
    lili=0
    #for i in range(1,2*n):
    while alin[-1]>0 and step<1.3*n:
        step=step+1
        ran=random.randint(1,alin[-1])
        k=bisect.bisect_left(alin,ran)
        link1=random.choice(range(len(A[k])))
        link2=random.choice(range(link1+2,link1+2+len(A[k])-3))%len(A[k])
        #pp=[np.exp(-gamma*min(m,len(A[k])-2-m)) for m in range(1,len(A[k])-2)]
        #link2=np.random.choice(range(link1+2,link1+2+len(A[k])-3),p=np.divide(pp,sum(pp)))%len(A[k])
        encircling=False
        addlink(A[k][link1],A[k][link2],A,H,tried)
        alen=[len(ii)*(len(ii)-3)/2 for ii in A]
        alin=[sum(alen[0:ii]) for ii in range(1,len(A)+1)]
        Als=sum([len(x)*(len(x)-3)/2 for x in A])
        ll=sum([(-1.5+np.sqrt(1.5**2+2*ii)+3)*alen.count(ii)*ii/float(Als) for ii in looplens])
        Links=list(set(H.edges()) - set(nx.cycle_graph(n).edges()))
    #print 'total time: ',(time.time()-starttime)/60.0,'min'
    #print 'alinmax',alen
    #return [n,nx.diameter(H)]
    #return np.sqrt(max(alin))
    #return nx.diameter(H)
    return Links
#print dias
#np.savetxt('takethosen_large',np.logspace(3.5,4, num=75))
#pp=[np.exp(-1*min(m,100-2-m)) for m in range(1,100-2)]
#print np.divide(pp,sum(pp))

#H=Fold_my_chain(500,0.1)
#print H.edges()
rest=[]
restgamma=[]
n=1000
Li=Fold_my_chain(n,0.00001)
#Len=[Li[i][1]-Li[i][0] for i in range(len(Li))]
Len=[min(abs(Li[i][1]-Li[i][0]),abs(Li[i][1]-Li[i][0]-n)) for i in range(len(Li))]
#print Len
Luna=[Len.count(i)/(1.0) for i in range(2,n/2)]
#print Luna
np.savetxt('luna_noencircle10',Luna)

print 'Luna sum: ',sum(Luna)

plotx=np.logspace(0.3, 1.7, num=12)
#bin_means, bin_edges, binnumber = sp.stats.binned_statistic(np.array(range(n/2-2)),np.array(Luna), statistic=np.mean, bins=plotx)
#print bin_means
#print bin_edges
##fit1, fitcov1=curve_fit(func, range(2,n/2),Luna)
##fit1, fitcov1=curve_fit(func, bin_edges[1:20],bin_means[0:20])
#fit1, fitcov1=curve_fit(func, bin_edges[1:20], bin_means)
bin_means, bin_edges, binnumber = sp.stats.binned_statistic(np.array(range(n/2-2)),np.array(Luna), statistic=np.mean, bins=plotx)
bin_centers = bin_edges[1:] - (bin_edges[1:] - bin_edges[0:11])/2
print bin_means
print bin_centers
#fit1, fitcov1=curve_fit(func, range(2,n/2),Luna)
#fit1, fitcov1=curve_fit(func, bin_edges[1:20],bin_means[0:20])
fit1, fitcov1=curve_fit(func, bin_centers[:], bin_means[:])

print fit1
#print np.logspace(0.01, 2.5, 50)

#lili=0
#for gamma in np.logspace(-2,1, num=20):
    #r1=[]
    #timer = time.time()
    #for lol in range(20):
        #r1.append(Fold_my_chain(1000,gamma))
    #rest.append(np.mean(r1))
    #restgamma.append(gamma)
    #print gamma,np.mean(r1),time.time()-timer
#resti1=[restgamma,rest]
#np.savetxt('gamma1000.dat', resti1)

#rest3=[]
#restgamma3=[]
#for gamma in np.logspace(-4,1, num=20):
    #r1=[]
    #timer = time.time()
    #for lol in range(20):
        #r1.append(Fold_my_chain(1000,gamma))
    #rest3.append(np.mean(r1))
    #restgamma3.append(gamma)
    #print gamma,np.mean(r1),time.time()-timer
#resti3=[restgamma3,rest3]
#np.savetxt('gamma_dia_1000_4.dat', resti3)

#resti1=np.loadtxt('gamma_dia_200_4.dat')
#resti2=np.loadtxt('gamma_dia_500_4.dat')
#resti3=np.loadtxt('gamma_dia_1000_4.dat')

#Plot linklength distribution
plt.plot(Luna)

plt.plot(range(2,n/2),func(np.array(range(2,n/2)),fit1[0], fit1[1]),label='fit: a=%.2f, b=%.2f'%(fit1[0],fit1[1]),color='red',linewidth=1)
plt.errorbar(bin_edges[1:],bin_means, marker='o',linewidth=0,label='data',color='black',elinewidth=1,barsabove=True)
plt.xlabel('$Linklength$', fontsize=20)
plt.ylabel('$Frequency$', fontsize=20)

#Plot weird phase transition things with gamma
#plt.plot(resti1[0],np.divide(resti1[1],1),label='chainlen 100',marker='*',linewidth=1, color='red')
#plt.plot(resti2[0],np.divide(resti2[1],1),label='chainlen 200',marker='*',linewidth=1, color='green')
#plt.plot(resti3[0],np.divide(resti3[1],1),label='chainlen 500',marker='*',linewidth=1, color='blue')
#plt.legend(loc=2,prop={'size':15})
#plt.xlabel('$\gamma$', fontsize=20)
#plt.ylabel('$R_g$', fontsize=20)
#nx.draw(H,node_size=10)
#plt.plot(dias.T[0],dias.T[1],label='simulation',marker='*',linewidth=0)

plt.yscale('log')
plt.xscale('log')
plt.show()