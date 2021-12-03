# -*- coding: utf-8 -*-
from __future__ import print_function

import glob
import os

import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import scipy as sp
import scipy.stats
import seaborn as sns
from Bio.PDB import *
from MDAnalysis import *
from scipy.optimize import curve_fit
import pickle

sns.set(font_scale=1.5)  # crazy big
sns.set_style("ticks")
sns.despine()


def func(x, a, b):
    return a * x ** (-b)


def func2(x, a):
    return a / x


def theofunc(x, C, N, ex):
    return (C) * np.power(N / (2 * (-x ** 2 + N * x - N)), ex)


def connected_component_subgraphs(G):
    for c in nx.connected_components(G):
        yield G.subgraph(c)


# def theofunc(x, C,N):
#     return (C)* N/(2*(N*x-N))

names_xray = glob.glob("../data/pdb_data2/all_xray/*.pdb")
names_xray = glob.glob("../data/pdb_data2/xray_data_short/*.pdb")
names_xray = glob.glob("../data/pdb_data2/xray_data/*.pdb")

names_nmr = glob.glob("../datapdb_data2/nmr_data/*.pdb")
names_nmr = glob.glob("../data/pdb_data2/nmr_data_short/*.pdb")

names_both = glob.glob("../data/pdb_data2/both/*.pdb")

# names_both = glob.glob("idp/*.pdb")

names_used = names_both

num = 5
lengths = [3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 16, 17, 19, 21, 23, 25, 28, 30, 33, 36, 40, 44, 48, 52, 58, 63, 69,
           76, 83, 91, 100, 110, 120, 132, 145, 158, 174, 191, 209, 229, 251, 275,
           302]#, 331, 363, 398]  # [0:45]
print(len(lengths))

numnum = 20
logscale = np.logspace(np.log10(3), np.log10(1000), num=numnum)
print('logscale', logscale)
results = []
binnedresults = []

if __name__ == '__main__':
    directories = ['../data/pdb_data2/both']
    # directories = ['idp']
    nameslen = len(names_used)
    # nameslen=50
    path = os.getcwd()
    meanies = []
    linklengths = []
    for d in directories:
        print(d)
        output = os.path.join(path, d)
        count = 0
        n = 1000
        t = 8.0
        counter = 0
        for counter in range(nameslen):
            # pdb_file = os.path.join(output,pdb.strip())+'.pdb'
            pdb_file = names_used[counter]
            print(pdb_file)
            u = Universe(pdb_file)
            # calphas = u.select_atoms("name CA and segid A")
            calphas = u.select_atoms("name CA and segid " + u.residues.segments.segids[0])
            print(counter, 'pdb file %s has %d calpha atoms' % (pdb_file, calphas.n_atoms))
            # coordinates are then easily accessed as
            r = calphas.positions
            maxdist = 0
            G = nx.empty_graph(len(r))
            avd = 0
            for i in range(len(r)):
                for j in range(i):
                    dist = np.linalg.norm(r[i] - r[j])
                    maxdist = max(dist, maxdist)
                    avd += dist ** 2
                    if dist < t:
                        G.add_edge(i, j)
                        linklengths.append(dist)
            avd = np.sqrt(avd) / (len(r))
            # take largest component
            if len(G.edges()) > 0:
                G = list(connected_component_subgraphs(G))[0]
                lena = len(G.nodes())
                Li = list(set(G.edges()) - set(nx.cycle_graph(lena).edges()))
                print(lena, 'nodes', len(Li), 'links')
                lililist = []
                for link in Li:
                    lililist.append(abs(link[0] - link[1]))
                results.append([lena, np.histogram(lililist, bins=range(3, 1000))])

print(count, nameslen)
# np.savetxt('/home/nora/Desktop/research/Göttingen/protein folding/3d/meaniesnewidp.txt',meanies)
# Saving pdb analysis
with open('pdb_results.pkl', 'wb') as f:
    pickle.dump(results, f)

##read in graph from adjacency matrix
sim_results = []
for k in lengths:
    length_list_sim = []
    lililist = []  # making a histogram of all sims of one length
    for i in np.random.permutation(range(1, num + 1)):
        M = np.array(
            np.loadtxt('../data/K50/Nora50_matrix_%d_%d.txt' % (k, i)))
        N = len(M)
        M[M > 1] = 0
        G = nx.from_numpy_matrix(M)
        sim_Li = list(G.edges())
        # print(k,'nodes',len(sim_Li),'links')
        # lililist=[] #making histogram of just one
        for link in sim_Li:
            lililist.append(abs(link[0] - link[1]))
    sim_results.append([k, np.histogram(lililist, bins=range(3, 1000))])

sns.set(font_scale=2)  # crazy big
sns.set_style("ticks")
sns.despine()
sns.set_palette(sns.xkcd_palette(["tomato", "windows blue", "amber", "greyish", "dusty purple", "faded green"]))

print('mean linklength', np.mean(linklengths))
heixdip = 2
countbin = 0
sum = np.zeros(996)
bsum = np.zeros(numnum - 1)
for i in range(len(results)):
    me = results[i]
    bin_means, bin_edges, binnumber = sp.stats.binned_statistic(me[1][1][:-1], me[1][0], statistic=np.mean,
                                                                bins=logscale)
    bin_centers = bin_edges[1:] - (bin_edges[1:] - bin_edges[0:numnum - 1]) / 2
    if me[0] < 260 and me[0] > 130:
        sum = sum + me[1][0]
        bsum = bsum + bin_means
        countbin += 1
print('number of len around 200 pdb', countbin)
print(me[1][0][0])
plt.plot(me[1][1][:-1], sum / countbin, linewidth=0, marker='o', label='PDB 200')
fit, ficov = curve_fit(func, me[1][1][heixdip:-1], sum[heixdip:] / countbin)
print('fit', fit)
bfit, bficov = curve_fit(func, bin_centers, bsum / countbin)
print('binned fit', bfit)
plt.plot(me[1][1][:-1], func(me[1][1][:-1], fit[0], fit[1]))

sum = np.zeros(996)
bsum = np.zeros(numnum - 1)
countbin = 0
for i in range(len(sim_results)):
    me = sim_results[i]
    bin_means, bin_edges, binnumber = sp.stats.binned_statistic(me[1][1][:-1], me[1][0], statistic=np.mean,
                                                                bins=logscale)
    bin_centers = bin_edges[1:] - (bin_edges[1:] - bin_edges[0:numnum - 1]) / 2
    if me[0] < 210 and me[0] > 200:
        sum = sum + me[1][0]
        bsum = bsum + bin_means
        countbin += 1
print('number of len around 200 sim', (countbin * num))
print(me[1][0][0])

simlen = np.linspace(0, 100, num=50)
print(simlen)
plt.plot(me[1][1][:-1], sum / (countbin * num), linewidth=0, marker='o', label='sim 200')
fit, ficov = curve_fit(func, me[1][1][heixdip:-1], sum[heixdip:] / (countbin * num))
print('sim fit', fit)
bfit, bficov = curve_fit(func, bin_centers, bsum / countbin)
print('sim binned fit', bfit)
plt.plot(me[1][1][:-1], func(me[1][1][:-1], fit[0], fit[1]))
plt.plot(np.array(simlen), theofunc(np.array(simlen), 80.0, 200.0, 0.8), label='theory')
plt.yscale('log')
plt.xscale('log')
plt.xlim([4, 100])
plt.ylim([0.01, 200])
plt.legend()
plt.savefig('both_200.pdf')

plt.figure()
countbin = 0
sum = np.zeros(996)
bsum = np.zeros(numnum - 1)
for i in range(len(results)):
    me = results[i]
    bin_means, bin_edges, binnumber = sp.stats.binned_statistic(me[1][1][:-1], me[1][0], statistic=np.mean,
                                                                bins=logscale)
    bin_centers = bin_edges[1:] - (bin_edges[1:] - bin_edges[0:numnum - 1]) / 2
    if me[0] < 460 and me[0] > 330:
        # plt.plot(me[1][1][:-1],me[1][0],label='%.0f'%me[0])
        # plt.plot(me[1][1][:-1],me[1][0],linewidth=0, marker='o', color='grey')
        sum = sum + me[1][0]
        bsum = bsum + bin_means
        countbin += 1
print('number of len around 400 pdb', countbin)
print(me[1][0][0])
plt.plot(me[1][1][:-1], sum / countbin, linewidth=0, marker='o', label='PDB 400')
# plt.plot(bin_centers,bsum/countbin, linewidth=0, marker='o')
fit, ficov = curve_fit(func, me[1][1][heixdip:-1], sum[heixdip:] / countbin)
print('fit', fit)
bfit, bficov = curve_fit(func, bin_centers, bsum / countbin)
print('binned fit', bfit)
plt.plot(me[1][1][:-1], func(me[1][1][:-1], fit[0], fit[1]))
sum = np.zeros(996)
bsum = np.zeros(numnum - 1)
countbin = 0
for i in range(len(sim_results)):
    me = sim_results[i]
    bin_means, bin_edges, binnumber = sp.stats.binned_statistic(me[1][1][:-1], me[1][0], statistic=np.mean,
                                                                bins=logscale)
    bin_centers = bin_edges[1:] - (bin_edges[1:] - bin_edges[0:numnum - 1]) / 2
    if me[0] < 430 and me[0] > 370:
        sum = sum + me[1][0]
        bsum = bsum + bin_means
        countbin += 1
print('number of len around 400 sim', (countbin * num))
print(me[1][0][0])

simlen = np.linspace(0, 200, num=100)
print(simlen)
plt.plot(me[1][1][:-1], sum / (countbin * num), linewidth=0, marker='o', label='sim 400')
fit, ficov = curve_fit(func, me[1][1][heixdip:-1], sum[heixdip:] / (countbin * num))
print('sim fit', fit)
bfit, bficov = curve_fit(func, bin_centers, bsum / countbin)
print('sim binned fit', bfit)
plt.plot(me[1][1][:-1], func(me[1][1][:-1], fit[0], fit[1]))
plt.plot(np.array(simlen), theofunc(np.array(simlen), 80.0, 400.0, 0.7), label='theory')
plt.yscale('log')
plt.xscale('log')
plt.xlim([4, 200])
plt.ylim([0.01, 200])
plt.legend()
plt.savefig('both_400.pdf')

plt.show()
