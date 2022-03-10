import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

alphafold_data = pd.read_csv("/Users/jasminguven/Documents/GitHub/sequence_distance_distribution/data/pdb_statistics"
                             "/AlphaFold_used.csv")
alpha_bins = alphafold_data["Bins"].to_numpy()
alpha_frequency = alphafold_data["Number"].to_numpy()

swissprot_data = pd.read_csv("/Users/jasminguven/Documents/GitHub/sequence_distance_distribution/data/pdb_statistics"
                             "/UniProt_Swiss_Prot.csv")
swiss_bins = swissprot_data["Bins"].to_numpy()
swiss_frequency = swissprot_data["Number"].to_numpy()
swiss_adjusted_frequency = swiss_frequency - alpha_frequency

# trembl_data = pd.read_csv("/Users/jasminguven/Documents/GitHub/sequence_distance_distribution/data/pdb_statistics"
#                           "/UniProt_TrEMBL_sorted.csv")
# trembl_bins = trembl_data["Bins"].to_numpy()
# trembl_frequency = trembl_data["Number"].to_numpy()
# trembl_adjusted_frequency = trembl_frequency - swiss_frequency - alpha_frequency

ll_data = pd.read_csv("/Users/jasminguven/Documents/GitHub/sequence_distance_distribution/data/pdb_statistics"
                      "/ll_used.csv")
ll_bins = ll_data["Bins"].to_numpy()
ll_frequency = ll_data["Number"].to_numpy()

rcsb_data = pd.read_csv("/Users/jasminguven/Documents/GitHub/sequence_distance_distribution/data/pdb_statistics"
                        "/RCSB_by_length.csv")

rcsb_bins = rcsb_data["Bins"].to_numpy()
rcsb_frequency = rcsb_data["Number"].to_numpy() - ll_frequency

sns.set_theme(context="notebook", palette="colorblind", style="ticks", font_scale=1.5, font="Helvetica")
fig = plt.figure(figsize=(10, 8))
ax = fig.subplots(2, 1, sharex=True)
# ax[0].bar(rcsb_bins, trembl_adjusted_frequency, color="#006374", label="TrEMBL PDBs", alpha=0.8, bottom=swiss_frequency)
ax[0].bar(rcsb_bins, swiss_adjusted_frequency, color="#006374", label="Swiss-Prot sequences",
          bottom=alpha_frequency)
# ax[0].bar(rcsb_bins, alpha_frequency, color="#fbafe4", label=r"Used alphaFold sequences")

ax[0].tick_params(axis="x", labelrotation=90)
# ax[0].set_yscale("log")
ax[0].legend()

ax[1].bar(rcsb_bins, rcsb_frequency, color="#006374", label="RCSB sequences", bottom=ll_frequency)
ax[1].bar(rcsb_bins, ll_frequency, color="#fbafe4", label="Used RCSB sequences")
ax[1].tick_params(axis="x", labelrotation=90)

fig.text(0.5, 0.025, "Chain length", ha="center")
plt.subplots_adjust(left=0.09, bottom=0.08, top=0.99, wspace=0.05, right=1)

ax[1].legend()
plt.tight_layout()
sns.despine()
plt.show()
