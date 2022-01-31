import seaborn as sbn
import networkx as nx
import pickle
import scipy as sp
import scipy.stats
from scipy.optimize import curve_fit
import glob
import os
import Bio
import MDAnalysis as mda
import numpy as np 
import warnings
warnings.filterwarnings('ignore')
import matplotlib.pyplot as plt
from sklearn.metrics import r2_score
import scipy.stats as stats
import pylab
def nth_harmonic(n):
    '''
    Compute the nth harmonic number, i.e. the sum of the reciprocals of the first n natural numbers.

    Parameters
    ----------
    n: int
        The number of natural numbers for the harmonic.

    Return
    ------
    harmonic: float
        The nth harmonic number.
    '''
    # Ensure n is an integer
    N = int(n)
    # Initialise the harmonic number as 1
    harmonic = 1.00
    for i in range(2, N + 1):
        harmonic += 1/i
    return harmonic


def f_high(link_length, chain_length, half_N_harmonic):
    '''
    Calculate f for k >> 1.

    Parameters
    ----------
    link_length: int
        The separation s betweeen each link. Bounds: 2 <= s < N/2
    chain_length: int
        The chain length, i.e. number N of C-alphas in a chain.
        Also the maxium possible links at the start.
    half_N_harmonic: float
        The "N/2"-th harmonic number.

    Return
    ------
    f_high: float
        The sequence distribution evaluated at a given link_length,
        with k >> 1.
    '''
    H_s = nth_harmonic(link_length)
    #f_high = (2/chain_length) * (((half_N_harmonic - H_s + 1) / (half_N_harmonic - 1)) * link_length - (half_N_harmonic / (half_N_harmonic - 1)))
    f_high = 2/chain_length *((half_N_harmonic - H_s + 1)/(half_N_harmonic-1)*link_length- half_N_harmonic/(half_N_harmonic-1))
    return f_high


def f_low(link_length, chain_length):
    '''
    Calculate f for k << N/2.

    Parameters
    ----------
    link_length: int
        The separation s betweeen each link. Bounds: 2 <= s < N/2
    chain_length: int
        The chain length, i.e. number N of C-alphas in a chain.
        Also the maxium possible links at the start.

    Return
    ------
    f_low: float
        The sequence distribution evaluated at a given link_length,
        with k << N/2.
    '''
    #f_low = - (2/chain_length**2) * link_length**2 + ((2/chain_length) + (6/chain_length))*link_length - ((2/chain_length) + (4/chain_length**2))
    f_low = -2/chain_length**2*link_length**2+(2/chain_length+6/chain_length**2)*link_length-(2/chain_length+4/chain_length**2)
    return f_low

def P_link_lengths(link_length, chain_length, half_N_harmonic, a, A):
    '''
    Evaluate the probability distribution of the realised residue distances
    of all added links, as given by the final approximation in Eq.(12) with
    all C_k = A.

    Parameters
    ----------
    link_length: int
        The separation s betweeen each link. Bounds: 2 <= s < N/2
    chain_length: int
        The chain length, i.e. number N of C-alphas in a chain.
        Also the maxium possible links at the start.
    half_N_harmonic: float
        The "N/2"-th harmonic number.
    a: int
        The number of steps used from Eq.(7).
    A: int
        The factor from the geometric series in Eq.(12), i.e. C_k = A.

    Return
    ------
    P: float
        Probability distribution of the realised residue distances 
        of all added links.
    '''
    f_l = f_low(link_length, chain_length)
    f_h = f_high(link_length, chain_length, half_N_harmonic)
    P = (((1-f_l)/(1-f_h))**a) * (A/f_h)
    return P

def P_in_link_lengths(link_length, chain_length, half_N_harmonic, a, k):
    '''
    Evaluate the probability distribution of the realised residue distances
    of all added links, as given by Eq.(12), with N/(N-3) approx. 1 and 
    C_k not all set to C.
    
    Parameters
    ----------
    link_length: int
        The separation s betweeen each link. Bounds: 2 <= s < N/2
    chain_length: int
        The chain length, i.e. number N of C-alphas in a chain. 
        Also the maxium possible links at the start. 
    half_N_harmonic: float
        The "N/2"-th harmonic number. 
    a: int
        The number of steps used from Eq.(7).
    k: int 
        The power k used in Eq.(11), i.e. the step number.
    
    Return
    ------
    P_in: float
        Probability distribution of the realised residue distances 
        of all added links.
    '''
    
    sum_range_C_k = range(2, N//2)
    
    for s in sum_range_C_k:
        f_low = f_low(s, chain_length)
        f_high = f_high(s, chain_length, half_N_harmonic)
        C_k = sum(F_k(f_low, f_high, chain_length, a, k))
    
    f_low = f_low(link_length, chain_length)
    f_high = f_high(link_length, chain_length, half_N_harmonic)
    P_in = (F_k(f_low, f_high, chain_length, a, k))/C_k 
    return P_in

def F_k(f_low, f_high, chain_length, a, k):
    '''
    Calculate the F_k distribution as given in Eq.(11). 
    That is, create the pool from which we draw links at step k. 
    
    Parameters
    ----------
    f_low: float
        The sequence distribution evaluated at a given link_length, 
        with k << N/2.
    f_high: float
        The sequence distribution evaluated at a given link_length, 
        with k >> 1.
    chain_length: int
        The chain length, i.e. number N of C-alphas in a chain. 
        Also the maxium possible links at the start. 
    a: int
        The number of steps used from Eq.(7).
    k: int
        The power k used in Eq.(11), i.e. the step number.

    Return
    ------
    F_k: float
        The function given in Eq.(11).
    '''
    F_k = chain_length * (((1 - f_low)/(1 - f_high))**a) * (1 - f_high)**k
    return F_k

def powerlaw(scale_factor, link_length_array, chain_length):
    '''
    Compute the power law relationship: N * (1/x^1.2)
    
    Parameters
    ----------
    link_length_array: array
        The separations s betweeen each link. Bounds: 2 <= s < N/2
    chain_length: int
        The chain length, i.e. number N of C-alphas in a chain. 
        Also the maxium possible links at the start.         
    
    Return
    ------
    pwr: float
        Power law relationship N * (1/x^1.2)
    '''
    pwr =  N * (1/pow(link_length_array * scale_factor, 1.2))
    return pwr


def get_simulation_data(length):
    '''
    Get simulation data for a given length.
    
    Parameters
    ----------
    length: str
        Link length: 100, 200 or 300.
    
    Return
    ------
    sim_data_tuple: tuple
        A tuple containing the simulation data histogram and the
        data edges from the histogram.
    '''
    simulation_data_files = glob.glob(f'../data/simulation_data/Nora50_matrix_{length}_*')
    link_length_list = []
    histogram = None
    # Loop over simulation files
    for file in simulation_data_files:
        matrix = np.loadtxt(file)
        '-----------------------'
        matrix[matrix > 1] = 0 # ???
        '-----------------------'
        protein_graph = nx.from_numpy_matrix(matrix)
        links = list(protein_graph.edges())
        for link in links:
            link_length_list.append(abs(link[0] - link[1]))
        histogram, edges = np.histogram(link_length_list, bins = range(3, 1000), density = True)
    # Normalise over the number of simulations
    histogram = histogram / len(simulation_data_files)
    sim_data_tuple = (histogram, edges)
    return sim_data_tuple

    
    
def normalise_data(data_array, bin_width):
    '''
    Take a data array and width of bins and normalise the data array.
    
    Parameters
    ----------
    data_array: np.ndarray
        Array of data to be normalised
    bin_width: int
        Width of a bin from the histogram.
   
    Return
    ------
    normalised_array: np.ndarray
        Normalised data array.
    '''
    data_sum = 0
    for i in range(len(data_array)):
        data_sum += data_array[i] * bin_width
    normalised_array = data_array / data_sum
    return normalised_array


link_length_ranges = ['100', '200', '300']
constant_A = [0.005, 0.001, 0.001] # check these
constant_a = [3, 3, 3] # check these

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

bs_stats = pd.read_csv('../data/rcsb_data/bootstrapped_100s_with_stats.csv')
# n_pdbs = len(np.load('../data/pdb_data/ids_100.npy'))
# bs_m = bs.melt().astype(np.float64)


# bs_stats = bs_m.groupby('variable', as_index = False).agg(mean = ('value', np.mean),
#                                                           lower_bound = ('value', lambda val: np.quantile(val, q = 0.05)),
#                                                           upper_bound = ('value', lambda val: np.quantile(val, q = 0.95)))

means = bs_stats['mean'].to_numpy()
lower_b = bs_stats['lower_bound'].to_numpy()
upper_b = bs_stats['upper_bound'].to_numpy()
var = bs_stats['variable'].to_numpy()

#print(var)
# sim_data, sim_edges = get_simulation_data('302')
# sim_bin_centres = sim_edges[:-1] + np.diff(sim_edges)/2
# sim_bin_width = sim_edges[1] - sim_edges[0]
# sim_data = normalise_data(sim_data, sim_bin_width)

normed_means = normalise_data(means, 1)
normed_lower_b = normalise_data(lower_b, 1)
normed_upper_b = normalise_data(upper_b, 1)

#plt.scatter(var,normed_means, s = 10, c = 'k', label = 'PDB 100s')

#plt.scatter(sim_edges, sim_data, s=10, c = 'r', label = 'SIM 100s')

# N = int(var[-1]+1)
N = int(var[-1]+1)
#print(N)
H_N_2 = nth_harmonic(N)


A = np.arange(0.0015, 0.0025, 0.0005)
a = np.arange(0, 4, 1)
fig = plt.figure()
ax = fig.subplots(len(a), len(A), sharex = True, sharey = True)
print(len(a), len(A))
for i in range(len(a)):
    for j in range(len(A)):

        ax[i][j].scatter(var, normed_means, s=10, c="k", label="AF PDB 300s")
        ax[i][j].fill_between(var, normed_lower_b, normed_upper_b, color="gray", alpha=0.6, label="95% C.L.")
        ax[i][j].plot(sumrange, [P_link_lengths(s, N, H_N_2, a[i], A[i]) for s in sumrange], c="#984ea3", label="Theory")
        ax[i][j].set_yscale("log")
        ax[i][j].set_xscale("log")
        ax[i][-1].set_ylabel(f'a = {a[i]}', fontsize = 13, rotation = 0, labelpad=21)
        ax[i][-1].yaxis.set_label_position('right')
        ax[0][j].set_title(f'A = {A[j]:2f}')

plt.legend()
fig.text(0.5, 0.025, 's / a.u. ', ha='center',fontsize=15.5) # shared x label
fig.text(0.005, 0.5, 'P(s)', va='center', rotation='vertical',fontsize=15.5) # shared y label
plt.subplots_adjust(left = 0.06, bottom = 0.08, top=0.95, wspace=0.1, right = 0.95)

fig = plt.figure()
ax = fig.subplots(len(a), len(A), sharex = True, sharey = True)

for row in range(len(a)):
    for col in range(len(A)):
        f = [P_link_lengths(s,N,H_N_2,a[row],A[col]) for s in sumrange]

        residuals = normed_means - f
        mean_of_residuals = np.mean(residuals)
        std_of_residuals = np.std(residuals)
        sum_of_residuals = np.sum(residuals)
        r_squared = r2_score(normed_means, f)
        #print(f'R^2 score: {r_squared}')
        #print(f'Residuals mean: {mean_of_residuals} and std: {std_of_residuals}, sum: {sum_of_residuals}')
        RSS = np.sum(residuals**2)
        ax[row][col].scatter(sumrange, residuals, s=10, marker = '.', c='red', label = 'Residuals')
        ax[row][col].hlines(0, sumrange[0], sumrange[-1], ls = '--', color='k', label = 'Zero')
        # plt.yscale('log')
        # plt.xscale('log')
        ax[row][col].plot([],[], ls=' ',label = f'RSS: {RSS:.4f}')
        ax[row][col].plot([],[], ls=' ',label = f'$R^2$: {r_squared:.4f}')
        ax[row][col].plot([],[], ls=' ',label = f'R mean: {mean_of_residuals:.4f}')
        # ax[row][col].set_ylim(-0.08,0.0)
        ax[row][-1].set_ylabel(f'a = {a[row]}', fontsize = 13, rotation = 0, labelpad=21)
        ax[row][-1].yaxis.set_label_position('right')
        ax[0][col].set_title(f'A = {A[col]:2f}')
        ax[row][col].legend()
fig.text(0.5, 0.025, 's / a.u. ', ha='center',fontsize=15.5) # shared x label
fig.text(0.005, 0.5, 'Residuals', va='center', rotation='vertical',fontsize=15.5) # shared y label
plt.subplots_adjust(left = 0.06, bottom = 0.08, top=0.95, wspace=0.1, right = 0.95)
#
plt.show()
