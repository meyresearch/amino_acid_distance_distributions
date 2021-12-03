

# Imports
import numpy as np
import networkx as nx
import traceback
import os
from Bio.PDB import *
from MDAnalysis import *
import warnings
import urllib
import pandas as pd
import socket
import random
import operator
warnings.filterwarnings('ignore')

def connected_component_subgraphs(G):
    for c in nx.connected_components(G):
        yield G.subgraph(c)

# Get data from the pdb
def do_analysis(calphas, t=8.0):
    link_length_list = []
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
    return link_length_list

def save_current_data(length, ll_list, id_list):
    if len(ll_list)>0:
        sim_data_hist, sim_data_edges = np.histogram(ll_list, bins=range(3, 1000), normed=True)
        sim_data_hist = sim_data_hist/len(id_list)

        fname1 = 'simd_data_hist_%s_curr.npy' % length
        fname2 = 'Ids_%s_curr.npy' % length
        np.save(fname2, id_list)
        np.save(fname1, sim_data_hist)

def save_bootstrapped_current_data(length,ll_list,id_list, n_samples = 1000):
    print (np.shape(ll_list), flush=True)
    if len(ll_list)>0:

        sim_data_hist_list = []
        sim_data_edges = None
        # Now we need to boostrap the list a bunch of times
        for n in range(n_samples):
            bts_ll_list = [random.choice(ll_list) for _ in ll_list]
            #bts_ll_list = np.array(bts_ll_list)
            bts_ll_list = [item for items in bts_ll_list for item in items]
            sim_data_hist, sim_data_edges = np.histogram(bts_ll_list, bins=range(3, 1000), normed=True)
            sim_data_hist = sim_data_hist/len(id_list)
            sim_data_hist_list.append(sim_data_hist)
        sim_data_hist_list = np.array(sim_data_hist_list)
        mean_list = sim_data_hist_list.mean(axis = 0)
        std_list = sim_data_hist_list.std(axis = 0)

        fname1 = 'sim_data_hist_mean_%s_curr.npy' % length
        fname2 = 'sim_data_hist_std_%s_curr.npy' % length
        fname3 = 'Ids_%s_curr.npy' % length
        np.save(fname1, mean_list)
        np.save(fname2, std_list)
        np.save(fname3, id_list)

def save_ll_lists(l_100, l_200, l_300, l_400, count):
    name1 = '100/ll_list_%s.npy' % count
    name2 = '200/ll_list_%s.npy' % count
    name3 = '300/ll_list_%s.npy' % count
    name4 = '400/ll_list_%s.npy' % count
    np.save(name1,l_100)
    np.save(name2,l_200)
    np.save(name3,l_300)
    np.save(name4,l_400)

if __name__ == "__main__":
    # Reading the pdb IDs
    socket.setdefaulttimeout(20)
    query = pd.read_csv('../data/combined_ids.csv')
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
    with open("log.txt", "w") as log:
        for pdb in query.columns:
            if counter % 5000 == 0:
                print("We are at entry %d/%d!" % (counter, n_pdbs), flush=True)
                save_bootstrapped_current_data('100', ll_100, hundred)
                save_bootstrapped_current_data('200', ll_200, twohundered)
                save_bootstrapped_current_data('300', ll_300, threehundred)
                save_bootstrapped_current_data('400', ll_400, fourhundred)
                save_ll_lists(ll_100,ll_200,ll_300,ll_400, counter)
            download = 'https://files.rcsb.org/download/%s.pdb' % pdb
            try:
                file_name = '../data/temp/'+pdb+'.pdb'
                #bio.retrieve_pdb_file(pdb,pdir='../data/tmep/', file_format='pbd')
                urllib.request.urlretrieve(download, file_name)
            except Exception:
                print("failed at %s" % pdb, flush=True)
                traceback.print_exc(file=log)
                continue 
            if os.path.isfile(file_name):
                u = Universe(file_name)
                calphas = u.select_atoms("name CA and segid " + u.residues.segments.segids[0])
                chain_len = len(calphas)
                #print(chain_len)
                if chain_len in range(85,115):
                    hundred.append(pdb)
                    ll_100.append(do_analysis(calphas))
                    #print(len(ll_100))
                if chain_len in range(185,215):
                    twohundered.append(pdb)
                    ll_200.append(do_analysis(calphas))
                if chain_len in range(285,315):
                    threehundred.append(pdb)
                    ll_300.append(do_analysis(calphas))
                if chain_len in range(385,415):
                    fourhundred.append(pdb)
                    ll_400.append(do_analysis(calphas))
                os.remove(file_name)
                counter +=1
            else:
                continue
    # data for length 100
    # sim_data_hist_100, sim_data_edges = np.histogram(ll_100, bins=range(3, 1000), normed=True)
    # sim_data_hist_100 = sim_data_hist_100/len(hundred)
    # np.save('simdata_hist_100.npy', sim_data_hist_100)
    # np.save('100_Ids.npy', hundred)

    # data for length 200
    # sim_data_hist_200, sim_data_edges = np.histogram(ll_200, bins=range(3, 1000), normed=True)
    # sim_data_hist_200 = sim_data_hist_200/len(twohundered)
    # np.save('simdata_hist_200.npy', sim_data_hist_200)
    # np.save('200_Ids.npy',twohundered)

    # data for length 300
    # sim_data_hist_300, sim_data_edges = np.histogram(ll_300, bins=range(3, 1000), normed=True)
    # sim_data_hist_300 = sim_data_hist_300/len(threehundred)
    # np.save('simdata_hist_300.npy', sim_data_hist_300)
    # np.save('300_Ids.npy',threehundred)

    # data for length 400
    # sim_data_hist_400, sim_data_edges = np.histogram(ll_400, bins=range(3, 1000), normed=True)
    # sim_data_hist_400 = sim_data_hist_400/len(fourhundred)
    # np.save('simdata_hist_400.npy', sim_data_hist_400)
    # np.save('400_Ids.npy',fourhundred)


