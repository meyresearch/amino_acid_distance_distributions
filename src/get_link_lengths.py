

# Imports
import numpy as np
import networkx as nx
import pickle
import scipy as sp
import scipy.stats
from scipy.optimize import curve_fit
import glob
import os
from Bio.PDB import *
from MDAnalysis import *
import warnings
import urllib
import pandas as pd
warnings.filterwarnings('ignore')



def connected_component_subgraphs(G):
    for c in nx.connected_components(G):
        yield G.subgraph(c)

# Get data from the pdb
def do_analysis(calphas, link_length_list, t=8.0):
    linklengths = []
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
        for link in Li:
            #return(abs(link[0]-link[1]))
            link_length_list.append(abs(link[0] - link[1]))

if __name__ == "__main__":
    # Reading the pdb IDs
    query = pd.read_csv('../data/mix')
    hundred = []
    twohundered = []
    threehundred = []
    fourhundred = []
    ll_100 = []
    ll_200 = []
    ll_300 = []
    ll_400 = []
    counter = 0
    n_pdbs = len(query.columns)
    for pdb in query.columns:
        if counter%1000 == 0:
            print("We are at entry %d/%d!" %(counter,n_pdbs))
        #print (pdb)
        download = 'https://files.rcsb.org/download/%s.pdb' %pdb
        #print(download)
        try:
            file_name = '../data/temp/'+pdb+'.pdb'
            urllib.request.urlretrieve(download, file_name)
        except:
            pass
        u = Universe(file_name)
        calphas = u.select_atoms("name CA and segid " + u.residues.segments.segids[0])
        chain_len = len(calphas)
        #print(chain_len)
        if chain_len in range(85,115):
            hundred.append(pdb)
            do_analysis(calphas, ll_100)
            #print(len(ll_100))
        if chain_len in range(185,215):
            twohundered.append(pdb)
            do_analysis(calphas, ll_200)
        if chain_len in range(285,315):
            threehundred.append(pdb)
            do_analysis(calphas, ll_300)
        if chain_len in range(385,415):
            fourhundred.append(pdb)
            do_analysis(calphas, ll_400)
        os.remove(file_name)
        counter +=1
    # data for length 100
    sim_data_hist_100, sim_data_edges = np.histogram(ll_100, bins=range(3, 1000), normed=True)
    sim_data_hist_100 = sim_data_hist_100/len(hundred)
    np.save('simdata_hist_100.npy', sim_data_hist_100)
    np.save('100_Ids.npy',threehundred)

    # data for length 200
    sim_data_hist_200, sim_data_edges = np.histogram(ll_200, bins=range(3, 1000), normed=True)
    sim_data_hist_200 = sim_data_hist_200/len(hundred)
    np.save('simdata_hist_200.npy', sim_data_hist_200)
    np.save('200_Ids.npy',twohundered)

    # data for length 300
    sim_data_hist_300, sim_data_edges = np.histogram(ll_300, bins=range(3, 1000), normed=True)
    sim_data_hist_300 = sim_data_hist_300/len(hundred)
    np.save('simdata_hist_300.npy', sim_data_hist_300)
    np.save('300_Ids.npy',threehundred)

    # data for length 400
    sim_data_hist_400, sim_data_edges = np.histogram(ll_400, bins=range(3, 1000), normed=True)
    sim_data_hist_400 = sim_data_hist_400/len(hundred)
    np.save('simdata_hist_400.npy', sim_data_hist_400)
    np.save('400_Ids.npy',fourhundred)


