import mdtraj as md
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import itertools

ref_structure_wt = "frame1_traiettoria_wt_paroxetina.gro"
trajectory_wt = "wt_paroxetina_it_traiettoria_center.xtc"

ref_structure_k = "frame1_k141e_paroxetina.gro"
trajectory_k = "k141e_paroxetina_it_traiettoria_center.xtc"

t_wt = md.load(trajectory_wt, top=ref_structure_wt)
t_k = md.load(trajectory_k, top=ref_structure_k)

paroxetine = [196]
residues = [i for i in range(0, 196)]

contact_pairs = list(itertools.product(paroxetine, residues))

# COMPUTE CONTACTS
contacts_paroxetine_wt, paroxetine_pairs = md.compute_contacts(t_wt, contact_pairs)
contacts_paroxetine_k, paroxetine_pairs = md.compute_contacts(t_k, contact_pairs)

# CONTACTS DATAFRAME
contacts_paroxetine_wt_df = pd.DataFrame(contacts_paroxetine_wt)
contacts_paroxetine_wt_df_median = contacts_paroxetine_wt_df.median(axis=0)

contacts_paroxetine_k_df = pd.DataFrame(contacts_paroxetine_k)
contacts_paroxetine_k_df_median = contacts_paroxetine_k_df.median(axis=0)

# MEDIAN
contacts_paroxetine_wt_df = pd.DataFrame(contacts_paroxetine_wt)
contacts_paroxetine_wt_df_median = contacts_paroxetine_wt_df.median(axis=0)

contacts_paroxetine_k_df = pd.DataFrame(contacts_paroxetine_k)
contacts_paroxetine_k_df_median = contacts_paroxetine_k_df.median(axis=0)

# DATAFRAME
paroxetine_distance_wt_df = pd.DataFrame({'Paroxetine': [i + 1 for i in paroxetine_pairs[:, 0]],
                                          'Residues': [i + 1 for i in paroxetine_pairs[:, 1]],
                                          'Distance': contacts_paroxetine_wt_df_median})

paroxetine_distance_k_df = pd.DataFrame({'Paroxetine': [i + 1 for i in paroxetine_pairs[:, 0]],
                                         'Residues': [i + 1 for i in paroxetine_pairs[:, 1]],
                                         'Distance': contacts_paroxetine_k_df_median})

# PIVOT
paroxetine_contacts_matrix_wt = paroxetine_distance_wt_df.pivot(index='Paroxetine', columns='Residues',
                                                                values='Distance')

paroxetine_contacts_matrix_k = paroxetine_distance_k_df.pivot(index='Paroxetine', columns='Residues',
                                                                values='Distance')

vmin = min(min(contacts_paroxetine_wt_df_median), min(contacts_paroxetine_k_df_median))
vmax = max(max(contacts_paroxetine_wt_df_median), max(contacts_paroxetine_k_df_median))

# HEATMAPS
fig, ax = plt.subplots(2, 1, sharex=True)

sns.heatmap(paroxetine_contacts_matrix_wt, cmap="Blues", ax=ax[0], vmin=vmin, vmax=vmax, cbar_kws={'label': 'Distance (nm)'})
sns.heatmap(paroxetine_contacts_matrix_k, cmap="Blues", ax=ax[1], vmin=vmin, vmax=vmax, cbar_kws={'label': 'Distance (nm)'})

ax[0].set_title("wt", fontsize=20)
ax[1].set_title("K141E", fontsize=20)

ax[0].set_xlabel(None)
ax[1].set_xlabel('Residue', fontsize=20)

ax[0].set_ylabel('Paroxetine', fontsize=20)
ax[1].set_ylabel('Paroxetine', fontsize=20)

plt.show()

fig.savefig("paroxetine_contactmaps_median_wt_k.png", dpi=320)








