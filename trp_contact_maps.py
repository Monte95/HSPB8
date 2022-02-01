import mdtraj as md
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import itertools
from trp import *

contact_pairs = list(itertools.product(trp, residues))

contacts_trp_wt, trp_pairs = md.compute_contacts(t_wt, contact_pairs)
contacts_trp_k, trp_pairs = md.compute_contacts(t_k, contact_pairs)

contacts_trp_wt_df = pd.DataFrame(contacts_trp_wt)
contacts_trp_wt_df_median = contacts_trp_wt_df.median(axis=0)

contacts_trp_k_df = pd.DataFrame(contacts_trp_k)
contacts_trp_k_df_median = contacts_trp_k_df.median(axis=0)

trp_distance_wt_df = pd.DataFrame({'TRP': [i+1 for i in trp_pairs[:, 0]],
                                   'Residues': [i+1 for i in trp_pairs[:, 1]],
                                   'Distance': contacts_trp_wt_df_median})

trp_distance_k_df = pd.DataFrame({'TRP': [i+1 for i in trp_pairs[:, 0]],
                                  'Residues': [i+1 for i in trp_pairs[:, 1]],
                                  'Distance': contacts_trp_k_df_median})

trp_contacts_matrix_wt = trp_distance_wt_df.pivot(index='TRP', columns='Residues', values='Distance')
trp_contacts_matrix_k = trp_distance_k_df.pivot(index='TRP', columns='Residues', values='Distance')

vmin = min(min(contacts_trp_wt_df_median), min(contacts_trp_k_df_median))
vmax = max(max(contacts_trp_wt_df_median), max(contacts_trp_k_df_median))


############
# HEATMAPS
############

fig, ax = plt.subplots(2, 1, sharex=True)

sns.heatmap(trp_contacts_matrix_wt, cmap="rocket_r", ax=ax[0], vmin=vmin, vmax=vmax, cbar_kws={'label': 'Distance (nm)'})
sns.heatmap(trp_contacts_matrix_k, cmap="rocket_r", ax=ax[1], vmin=vmin, vmax=vmax, cbar_kws={'label': 'Distance (nm)'})

ax[0].set_title("wt", fontsize=20)
ax[1].set_title("K141E", fontsize=20)

ax[0].set_xlabel(None)
ax[1].set_xlabel('Residue', fontsize=20)

ax[0].set_ylabel('Trp residue', fontsize=20)
ax[1].set_ylabel('Trp residue', fontsize=20)

plt.show()

fig.savefig("trp_contactmaps_median_wt_k.png", dpi=320)




