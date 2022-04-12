import numpy as np
import scipy.stats
import theory_functions
import matplotlib.pyplot as plt


lengths = ["100", "200", "300"]
for length in lengths:
    rcsb = np.load(f"../data/rcsb/histogram_{length}_not_normed.npy", allow_pickle=True)
    af = np.load(f"../data/alphafold/histogram_{length}_not_normed.npy", allow_pickle=True)

    mean_rcsb = theory_functions.get_measure_of_central_tendency(rcsb, "mean")
    mean_af = theory_functions.get_measure_of_central_tendency(af, "mean")

    normalised_rcsb_mean = mean_rcsb / np.sum(mean_rcsb)
    normalised_af_mean = mean_af / np.sum(mean_af)
    print(f"Length: {length}")
    ks = scipy.stats.ks_2samp(normalised_rcsb_mean, normalised_af_mean, alternative="two-sided", mode="auto")
    print(f"KS statistic {ks.statistic}, p-value {ks.pvalue}")
    ks_range = scipy.stats.ks_2samp(normalised_rcsb_mean[4:int(length)//2], normalised_af_mean[4:int(length)//2],
                                    alternative="two-sided", mode="auto")
    print(f"Range 4:{int(length)//2}, KS statistic {ks_range.statistic}, p-value {ks_range.pvalue}")

    x = np.linspace(1, 300, 300)[:-1]
    af_cdf = np.cumsum(normalised_af_mean)
    rcsb_cdf = np.cumsum(normalised_rcsb_mean)
    plt.plot(x, af_cdf, "r", "--")
    plt.plot(x, rcsb_cdf, "k")
    # plt.xlim(0, int(length)//2)
    plt.show()
