import networkx as nx
import numpy as np
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
     return a * x**b
def func2(x, a):
     #return a * x**meanp
     return a * x**0.60228
def pi(x):
    #return 1.0/3
    if (x<n/2.0): return float(4*x)/(3*(n+2*x)) +1.0/3
    else: return float(2)/3
    #return float(3*x+4+n)/(3*(n+x))
def pimin(x):
    return 1.0/3


def diameter_size(G):
    shortest_path_pairwise=nx.shortest_path(G)
    all_shortest_paths=[]

    for shortest_paths_from_src in shortest_path_pairwise.values():
        all_shortest_paths+=shortest_paths_from_src.values()

    diameter=max([len(path) for path in all_shortest_paths])
    paths_forming_diameter=[path for path in all_shortest_paths if len(path)==diameter]
    return len(set([node for path in paths_forming_diameter for node in path]))
#flattening the list and using set() to get rid of duplicate nodes that are in multiple diameters

print(  '---------------------------------------')
starttime = time.time()
s=1
maxi=30
#fitmin=8
#mini=8
#meandias=np.zeros((maxi-mini,2))
#meandias2=np.zeros((maxi-mini,2))
#meanvol=np.zeros((maxi-mini,2))
#meanp=np.zeros((maxi-mini,2))
meandias=np.zeros((maxi,2))
meandias2=np.zeros((maxi,2))
meanvol=np.zeros((maxi,2))
meanp=np.zeros((maxi,2))

#for nrange in range(mini,maxi):
nrange=0
print( np.logspace(1.5, 3.4, num=maxi))
for n2 in np.logspace(1.5,3.4, num=maxi):
    n=int(n2)
    #n=int(50+1.1**nrange)
    looplens=[i*(i-3)/2 for i in range(3,n+1)]
    #create cycle graph
    H=nx.cycle_graph(n)
    Diam=np.zeros((s))
    Diam2=np.zeros((s))
    Vol=np.zeros((s))
    p=np.zeros((n-2,s))

    #Di[nrange-mini,1]=1
    #Di[nrange-mini,0]=n
    #for ni in range(2,n-2):
        ##print Di[nrange-mini,1]
        #Di[nrange-mini,1]=Di[nrange-mini,1]+(Di[nrange-mini,1]+2*(ni+1-Di[nrange-mini,1])/(Di[nrange-mini,1])-3)/(ni+1)*pi(ni)

    for j in range(s):
        print( j)
        H=nx.cycle_graph(n)
        alin=[]
        A=nx.cycle_basis(H)
        for i in range(1,n-2):
            A=nx.cycle_basis(H)
            alen=[len(ii)*(len(ii)-3)/2 for ii in A]
            alin=[sum(alen[0:ii]) for ii in range(1,len(A)+1)]
            ran=random.randint(1,max(alin))
            k=bisect.bisect_left(alin,ran)
            link1=random.choice(range(len(A[k])))
            link2=random.choice(range(link1+2,link1+2+len(A[k])-3))%len(A[k])
            #if flatten(A).count(A[k][link1])<6 and flatten(A).count(A[k][link2])<6:
                #H.add_edge(A[k][link1],A[k][link2])
            H.add_edge(A[k][link1],A[k][link2])

            Als=sum([len(x)*(len(x)-3)/2 for x in A])
            #print alen
            ll=sum([(-1.5+np.sqrt(1.5**2+2*ii)+3)*alen.count(ii)*ii/float(Als) for ii in looplens])
            #print ll
            #ll=sum([len(x)*len(x)*(len(x)-3) for x in A])//Als
            p[i][j]=(4+ll)/float(3*ll)

        Diam[j]=nx.diameter(H)
        D=nx.Graph()
        if len(A)>1:
            D.add_nodes_from(range(len(A)))
            for l in range(len(A)):
                for i in range(l+1,len(A)):
                    #if len([val for val in A[l] if val in A[i]])>0:
                    if len([val for val in A[l] if val in A[i]])>1:
                        D.add_edge(l,i)
            Diam2[j]=nx.diameter(D)
            Vol[j]=diameter_size(D)

    #meandias[nrange-mini,0]=n
    #meandias[nrange-mini,1]=np.mean(Diam)
    #meandias2[nrange-mini,0]=n
    #meandias2[nrange-mini,1]=np.mean(Diam2)
    #meanvol[nrange-mini,0]=n
    #meanvol[nrange-mini,1]=np.mean(Vol)
    #meanp[nrange-mini,0]=n
    #meanp[nrange-mini,1]=np.mean(p)

    meandias[nrange,0]=n
    meandias[nrange,1]=np.mean(Diam)
    meandias2[nrange,0]=n
    meandias2[nrange,1]=np.mean(Diam2)
    meanvol[nrange,0]=n
    meanvol[nrange,1]=np.mean(Vol)
    meanp[nrange,0]=n
    meanp[nrange,1]=np.mean(p)
    print( 'nodes: ',str(n), 'time:',(time.time()-starttime)/60.0,' min')
    nrange=nrange+1

print( '------list of mean diameters--------')
print( 'this took:',(time.time()-starttime)/60, 'min')
meanp=np.mean(p)
print( 'mean p=',meanp)

#maxi=96
Di=np.zeros((maxi,2))
Dimin=np.zeros((maxi,2))
nrange=0
for n2 in np.logspace(1.5, 4, num=maxi):
    n=int(n2)
    Di[nrange,1]=1
    Di[nrange,0]=n
    Dimin[nrange,1]=1
    Dimin[nrange,0]=n
    pp=np.zeros((n-4,2))
    for ni in range(2,n-2):
        pp[ni-2,0]=ni
        pp[ni-2,1]=pi(ni)
        #print ni,'pi= ',pi(ni)
        Di[nrange,1]=Di[nrange,1]+(Di[nrange,1]+2*(ni+1-Di[nrange,1])/(Di[nrange,1])-4)/(ni+1)*pi(ni)
        Dimin[nrange,1]=Dimin[nrange,1]+(Dimin[nrange,1]+2*(ni+1-Dimin[nrange,1])/(Dimin[nrange,1])-4)/(ni+1)*pimin(ni)
    nrange=nrange+1

#Di=np.zeros((maxi-mini,2))
#for nrange in range(mini,maxi):
    #n=int(50+1.1**nrange)
    #Di[nrange-mini,1]=1
    #Di[nrange-mini,0]=n
    #for ni in range(2,n-2):
        ##print Di[nrange-mini,1]
        #Di[nrange-mini,1]=Di[nrange-mini,1]+(Di[nrange-mini,1]+2*(ni+1-Di[nrange-mini,1])/(Di[nrange-mini,1])-3)/(ni+1)*pi(ni)
        #
#print pp
plt.figure(figsize=(8,8))
plt.plot(np.mean(p,axis=1))
plt.plot(pp.T[0],pp.T[1])
plt.show()
#nlist=np.arange(int(len(p)))
#reslist=4.0/3.0*nlist/(len(p)+2.0*nlist)+1.0/3.0
#plt.plot(nlist[0:int(len(p)/2)+1],reslist[0:int(len(p)/2)+1],color='green')
##plt.xlim([0,len(p)/2.0])
#plt.plot(nlist[int(len(p)/2):len(p)],np.ones(len(p))[int(len(p)/2):len(p)]*2.0/3.0,color='green')
##plt.plot(nlist,np.ones(len(p))*0.602284,color='red')
##plt.xlim([len(p)/2.0,len(p)])
##plt.savefig('p.pdf')

#np.savetxt('saved.dat', meandias.T[1])
#fit1, fitcov1=curve_fit(func, meandias.T[0][fitmin:maxi],meandias.T[1][fitmin:maxi])
#fit2, fitcov2=curve_fit(func, meandias2.T[0][fitmin:maxi],meandias2.T[1][fitmin:maxi])
#fitv, fitcovv=curve_fit(func, meanvol.T[0][fitmin:maxi],meanvol.T[1][fitmin:maxi])
#fitana, ana=curve_fit(func2, meandias.T[0][fitmin:maxi],meandias.T[1][fitmin:maxi])
fit1, fitcov1=curve_fit(func, meandias.T[0],meandias.T[1])
fit2, fitcov2=curve_fit(func, meandias2.T[0],meandias2.T[1])
fitv, fitcovv=curve_fit(func, meanvol.T[0],meanvol.T[1])
fitana, ana=curve_fit(func2, meandias.T[0],meandias.T[1])
perr=np.sqrt(np.diag(fitcov1))

print( 'triangle diameter scaling: ',fit1)
print( perr)
print( 'loop tree scaling:         ',fit2)
print( 'Volume scaling:            ',fitv)
print( 'analytics fit:            ',fitana)

matplotlib.rc('xtick',labelsize=20)
matplotlib.rc('ytick',labelsize=20)

#plt.figure(figsize=(8,8))
##plt.plot(meandias.T[0],func2(meandias.T[0],fitana), label='analytic result')
#plot1=plt.plot(meandias.T[0],func(meandias.T[0],fit1[0], fit1[1]),label='fit')
#plot2=plt.plot(meandias.T[0],meandias.T[1],label='simulation',marker='*',linewidth=0)
##plt.plot(Di.T[0],Di.T[1],label='semi-analytic simulation')
##plt.plot(meandias.T[0],func(meandias.T[0],fit2[0], fit2[1]))
##plt.plot(meandias2.T[0],meandias2.T[1])
##plt.plot(meanvol.T[0],meanvol.T[1],label='Volume simulation')
##plt.plot(meandias.T[0],func(meandias.T[0],fitv[0], fitv[1]),label='Volume fit')
#plt.xlim([50+1.1**mini,50+1.1**maxi])
#plt.yscale('log')
#plt.xscale('log')
#plt.xlabel('N')
#plt.ylabel('D(N)')
#plt.legend(loc=4)
#plt.savefig('scaling_numerics.pdf')
print( meandias.T[0])
plotx=np.logspace(1.5, 4, num=maxi)
plt.figure(figsize=(10,8))
#plt.plot(meandias.T[0],func2(meandias.T[0],fitana), label='analytic result')
plot1=plt.plot(plotx,func(plotx,fit1[0], fit1[1]),label='fit: nu=%.2f'%fit1[1],color='red',linewidth=3)
plot1=plt.plot(plotx,func(plotx,fit1[0]-0.5*perr[0], 1/2.0),label='space filling',color='blue')
plot1=plt.plot(plotx,func(plotx,fit1[0]+0.5*perr[0], 3/4.0),label='self-avoiding walk',color='green')

plt.fill_between(plotx, func(plotx,fit1[0]-perr[0], fit1[1]-perr[1]), func(plotx,fit1[0]+perr[0], fit1[1]+perr[1]), facecolor='red', alpha=0.3)
#plot1=plt.plot(meandias.T[0],func(meandias.T[0],fit1[0]-perr[0], fit1[1]-perr[1]))
#plot1=plt.plot(meandias.T[0],func(meandias.T[0],fit1[0]+perr[0], fit1[1]+perr[1]))
plot2=plt.plot(meandias.T[0],meandias.T[1],label='simulation',marker='*',linewidth=0)
plt.plot(Di.T[0],Di.T[1],label='max-p simulation',color='black')
plt.plot(Dimin.T[0],Dimin.T[1],label='min-p simulation',color='grey')

#plt.plot(meandias.T[0],func(meandias.T[0],fit2[0], fit2[1]))
#plt.plot(meandias2.T[0],meandias2.T[1],label='simulation',marker='*',linewidth=0)
#plt.plot(meanvol.T[0],meanvol.T[1],label='Volume simulation')
#plt.plot(meandias.T[0],func(meandias.T[0],fitv[0], fitv[1]),label='Volume fit')
#plt.xlim([min(meandias.T[0]),max(meandias.T[0])])
plt.xlim([min(np.logspace(1.5, 4, num=maxi)),max(np.logspace(1.5, 4, num=maxi))])
plt.ylim([min(meandias.T[1]),200])

plt.yscale('log')
plt.xscale('log')
plt.xlabel('Chain length N', fontsize=25)
plt.ylabel('Diameter D(N)', fontsize=25)
plt.legend(loc=4,prop={'size':17})
plt.savefig('newsemi.pdf')
print( nx.diameter(D))

#plt.figure(figsize=(8,8))
#nx.draw_graphviz(D,node_size=10)
#plt.axis('off')
#plt.savefig('dual1.pdf')

#plt.figure(figsize=(8,8))
#nx.draw_graphviz(H,node_size=10)
#plt.axis('off')
#plt.savefig('folded1.pdf')

plt.show()

print( 'Done! enjoy!')
