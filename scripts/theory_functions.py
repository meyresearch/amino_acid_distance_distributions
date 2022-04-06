"""
Theoretical functions for plotting.
"""
import numpy as np
from sklearn.metrics import r2_score
import networkx as nx


def harmonic_number(n_numbers: int) -> float:
    """
    @param n_numbers: number of natural numbers
    @return: harmonic: nth harmonic number
    """
    n_int = int(n_numbers)
    harmonic = 1.00
    for i in range(2, n_int + 1):
        harmonic += 1 / i
    return harmonic


def f_high(amino_acid_distance: int, chain_length: int, half_n_harmonic_number: float) -> float:
    """
    Calculate f for k >> 1
    @param amino_acid_distance: separation between each amino acid, 2 <= s < N/2
    @param chain_length: number of C-alphas in a chain
    @param half_n_harmonic_number: the 'N/2'-th harmonic number
    @return: amino acid distribution for k >> 1
    """
    harmonic = harmonic_number(amino_acid_distance)
    constant_multiplier = 2 / chain_length
    s_term = ((half_n_harmonic_number - harmonic + 1) / (half_n_harmonic_number - 1)) * amino_acid_distance
    constant_term = half_n_harmonic_number / (half_n_harmonic_number - 1)
    return constant_multiplier * (s_term - constant_term)


def f_low(amino_acid_distance: int, chain_length: int) -> float:
    """
    Calculate f for k << N/2
    @param amino_acid_distance: separation between each link, 2 <= s < N/2
    @param chain_length: number of C-alphas in a chain
    @return: amino acid distance distribution for k << N/2
    """
    s_square_term = (-2 / chain_length ** 2) * amino_acid_distance ** 2
    s_term = (2 / chain_length + 6 / chain_length ** 2) * amino_acid_distance
    constant_term = (2 / chain_length) + (4 / chain_length ** 2)
    return s_square_term + s_term - constant_term


def amino_acid_distance_distribution(amino_acid_distances: np.ndarray, chain_length: int, half_n_harmonic_number: float,
                                     exponent: int, dimensionality: float) -> np.ndarray:
    """

    @param amino_acid_distances: separation between each amino acid 2 <= s < N/2
    @param chain_length: number of C-alphas
    @param half_n_harmonic_number: the 'N/2'-th harmonic number
    @param exponent: constant a
    @param dimensionality: dimensionality scaling constant (constant A)
    @return: probability distribution of the realised amino acid distances
    """
    distributions = []
    for s in amino_acid_distances:
        f_low_k = f_low(s, chain_length)
        f_high_k = f_high(s, chain_length, half_n_harmonic_number)
        if f_high_k == 0.0:
            distributions.append(0)
        else:
            distribution = (((1 - f_low_k) / (1 - f_high_k)) ** exponent) * (dimensionality / f_high_k)
            distributions.append(distribution)
    return np.asarray(distributions)


def power_law(distances: np.ndarray, power: float, constant_multiplier: float):
    """
    Plot a power law 1/s^a
    @param distances: amino acid distances
    @param power: power to raise to
    @param constant_multiplier: constant
    @return: np.ndarray
    """
    powerlaw = constant_multiplier * 1 / pow(distances, power)
    return np.asarray(powerlaw)


def plotting_statistics(dimensionality: float, exponent: int, n_points: int, half_n_harmonic_number: float,
                        plotting_sumrange: np.ndarray, normalised_measure: np.ndarray) -> tuple:
    """
    Calculate plotting statistics: residuals, mean of residuals, R^2 statistic and residual sum of squares
    @param dimensionality: dimensionality scaling constant (constant A)
    @param exponent: constant a
    @param n_points: number of datapoints
    @param half_n_harmonic_number: "N/2"-th harmonic number
    @param plotting_sumrange: range to sum over for plotting
    @param normalised_measure: normalised means of amino acid distance frequencies
    @return: tuple of different statistics
    """
    theory_function = amino_acid_distance_distribution(plotting_sumrange, n_points, half_n_harmonic_number,
                                                       exponent, dimensionality)
    residuals = normalised_measure - theory_function
    residuals_mean = np.mean(residuals)
    residuals_sum = np.sum(residuals)
    residual_sum_of_squares = residuals_sum ** 2
    r_square_value = r2_score(normalised_measure, theory_function)
    print("------------------------------------------")
    print(f"a: {exponent}, A: {dimensionality}")
    print(f"mean of residuals: {residuals_mean}")
    print(f"sum of residuals: {residuals_sum}")
    print(f"r-squared statistic: {r_square_value}")
    print(f"residual sum of squares: {residual_sum_of_squares}")
    return residuals, residuals_mean, r_square_value, residual_sum_of_squares


def get_adjacency_matrix(protein_graph: nx.Graph) -> np.ndarray:
    """
    Take a protein graph and convert it to an adjacency matrix (Numpy array)
    @param protein_graph: protein contact map (graph)
    @return: adjacency matrix
    """
    return nx.convert_matrix.to_numpy_array(protein_graph)
