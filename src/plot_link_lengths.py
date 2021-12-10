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

def powerlaw(link_length_array, chain_length):
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
    pwr = N * (1/pow(link_length_array, 1.2))
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
constant_A = [0.001, 0.001, 0.001] # check these
constant_a = [3, 3, 3] # check these
#shift_x = [9, 8.5, 8.5] 
#shift_y = [500, 1000, 1000]
for ll_index in range(len(link_length_ranges)):
    
    # Load PDB data
    pdb_ids = np.load(f'../data/data_for_plotting/ids_{link_length_ranges[ll_index]}.npy')
    pdb_mean_data = np.load(f'../data/data_for_plotting/means_{link_length_ranges[ll_index]}.npy')
    
    # Get PDB bin widths and centres
    pdb_edges = np.load(f'../data/data_for_plotting/edges_{link_length_ranges[ll_index]}.npy')[0]
    pdb_bin_centres = pdb_edges[:-1] + np.diff(pdb_edges)/2
    pdb_bin_width = pdb_edges[1] - pdb_edges[0]
    
    # Normalise PDB data
    pdb_mean_data = normalise_data(pdb_mean_data, pdb_bin_width)
    
    # Load simulation data
    if link_length_ranges[ll_index] == '200':
        sim_data, sim_edges = get_simulation_data('209')
    elif link_length_ranges[ll_index] == '300':
        sim_data, sim_edges = get_simulation_data('302')
    else: 
        sim_data, sim_edges = get_simulation_data(link_length_ranges[ll_index])

    # Get simulation bin widths and centres
    sim_bin_centres = sim_edges[:-1] + np.diff(sim_edges)/2
    sim_bin_width = sim_edges[1] - sim_edges[0]
    
    # Normalise simulation data
    sim_data = normalise_data(sim_data, sim_bin_width)
    
    # Plotting
    fig = plt.figure(figsize = (8,4))
    ax = fig.add_subplot()
    
    # Plot PDB data
    ax.scatter(pdb_bin_centres, pdb_mean_data, marker = '.', c='#e41a1c',label = f'PDB {link_length_ranges[ll_index]} means')

    # Plot simulation data 
    ax.scatter(sim_bin_centres, sim_data, marker = '^', label = f'SIM {link_length_ranges[ll_index]}', c = '#999999')
    
    # Plot theory 
    N=int(link_length_ranges[ll_index])
    H_N_2=nth_harmonic(N//2)
    A=constant_A[ll_index]
    a=constant_a[ll_index]
    sumrange=range(2,N//2)
    plotrange=range(20,N//2)
    ax.plot(np.array(sumrange),[P_link_lengths(s,N,H_N_2,a,A) for s in sumrange],c='#984ea3',label='Theory')
    ax.plot(np.array(sumrange), powerlaw(np.array(sumrange), N), c='black',label = 'Power law')
    
    # Plot 97% confidence interval      
    means_sorted = np.sort(pdb_mean_data)
    lower_bound = means_sorted[len(pdb_mean_data)//3]
    upper_bound = means_sorted[len(pdb_mean_data)*2//3]        
    ax.fill_between(pdb_bin_centres, (pdb_mean_data-lower_bound), (pdb_mean_data+upper_bound), color = 'red', alpha = 0.2, label = '3$\sigma$ C.I.', zorder = -1)

    ax.set_yscale('log')
    ax.set_xscale('log')
    ax.set_xlim([2,500])    
    ax.set_ylim([0, 10])
    
    ax.legend(loc='lower left')
    ax.set_ylabel('P(s)', fontsize = 15)
    ax.set_xlabel('s / a.u.', fontsize = 15)
    plt.show()
    #plt.savefig(f'../data/{link_length_ranges[ll_index]}_plot.svg', format = 'svg')




