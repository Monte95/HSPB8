import mdtraj as md
import numpy as np
import pandas as pd
import seaborn as sns
import math
import matplotlib.pyplot as plt
import itertools
from mpl_toolkits.mplot3d import Axes3D

ref_structure_all_wt = "it_wt.gro"
trajectory_all_wt = "all_wt.xtc"

ref_structure_it_wt = "it_wt.gro"
trajectory_it_wt = "it_wt.xtc"

ref_structure_m_wt = "m_wt.gro"
trajectory_m_wt = "m_wt.xtc"

ref_structure_r_wt = "trros_wt.gro"
trajectory_r_wt = "trros_wt.xtc"

trajectory_emap_wt = ['generated_encodermap_wt_7.pdb', 'generated_encodermap_wt_8.pdb', 'generated_encodermap_wt_9.pdb',
                      'generated_encodermap_wt_10.pdb', 'generated_encodermap_wt_11.pdb',
                      'generated_encodermap_wt_12.pdb',
                      'generated_encodermap_wt_13.pdb', 'generated_encodermap_wt_14.pdb',
                      'generated_encodermap_wt_15.pdb',
                      'generated_encodermap_wt_16.pdb', 'generated_encodermap_wt_17.pdb',
                      'generated_encodermap_wt_18.pdb',
                      'generated_encodermap_wt_19.pdb']

all_wt = md.load(trajectory_all_wt, top=ref_structure_all_wt)
it_wt = md.load(trajectory_it_wt, top=ref_structure_it_wt)
m_wt = md.load(trajectory_m_wt, top=ref_structure_m_wt)
r_wt = md.load(trajectory_r_wt, top=ref_structure_r_wt)
emap_wt = md.load(trajectory_emap_wt)

contact_pairs = [(9, 98), (9, 194), (98, 194)]

distance_cys_all_wt, cys_pairs = md.compute_contacts(all_wt, contact_pairs)
distance_cys_it_wt, cys_pairs = md.compute_contacts(it_wt, contact_pairs)
distance_cys_m_wt, cys_pairs = md.compute_contacts(m_wt, contact_pairs)
distance_cys_r_wt, cys_pairs = md.compute_contacts(r_wt, contact_pairs)
distance_cys_emap_wt, cys_pairs = md.compute_contacts(emap_wt, contact_pairs)

# CYS10-CYS99
distance_cys_10_99_all_wt_df = pd.DataFrame({'Distance': distance_cys_all_wt[:, 0],
                                             'Homology Model': np.repeat('All wt', len(distance_cys_all_wt[:, 0]))})

distance_cys_10_99_it_wt_df = pd.DataFrame({'Distance': distance_cys_it_wt[:, 0],
                                            'Homology Model': np.repeat('I-TASSER', len(distance_cys_it_wt[:, 0]))})

distance_cys_10_99_m_wt_df = pd.DataFrame({'Distance': distance_cys_m_wt[:, 0],
                                           'Homology Model': np.repeat('MODELLER', len(distance_cys_m_wt[:, 0]))})

distance_cys_10_99_r_wt_df = pd.DataFrame({'Distance': distance_cys_r_wt[:, 0],
                                           'Homology Model': np.repeat('ROSETTA', len(distance_cys_r_wt[:, 0]))})

distance_10_99_df = pd.concat([distance_cys_10_99_all_wt_df, distance_cys_10_99_it_wt_df, distance_cys_10_99_m_wt_df,
                               distance_cys_10_99_r_wt_df], ignore_index=True)

# CYS10-CYS195
distance_cys_10_195_all_wt_df = pd.DataFrame({'Distance': distance_cys_all_wt[:, 1],
                                              'Homology Model': np.repeat('All wt', len(distance_cys_all_wt[:, 1]))})

distance_cys_10_195_it_wt_df = pd.DataFrame({'Distance': distance_cys_it_wt[:, 1],
                                             'Homology Model': np.repeat('I-TASSER', len(distance_cys_it_wt[:, 1]))})

distance_cys_10_195_m_wt_df = pd.DataFrame({'Distance': distance_cys_m_wt[:, 1],
                                            'Homology Model': np.repeat('MODELLER', len(distance_cys_m_wt[:, 1]))})

distance_cys_10_195_r_wt_df = pd.DataFrame({'Distance': distance_cys_r_wt[:, 1],
                                            'Homology Model': np.repeat('ROSETTA', len(distance_cys_r_wt[:, 1]))})

distance_10_195_df = pd.concat(
    [distance_cys_10_195_all_wt_df, distance_cys_10_195_it_wt_df, distance_cys_10_195_m_wt_df,
     distance_cys_10_195_r_wt_df], ignore_index=True)

# CYS99-CYS195
distance_cys_99_195_all_wt_df = pd.DataFrame({'Distance': distance_cys_all_wt[:, 2],
                                              'Homology Model': np.repeat('All wt', len(distance_cys_all_wt[:, 2]))})

distance_cys_99_195_it_wt_df = pd.DataFrame({'Distance': distance_cys_it_wt[:, 2],
                                             'Homology Model': np.repeat('I-TASSER', len(distance_cys_it_wt[:, 2]))})

distance_cys_99_195_m_wt_df = pd.DataFrame({'Distance': distance_cys_m_wt[:, 2],
                                            'Homology Model': np.repeat('MODELLER', len(distance_cys_m_wt[:, 2]))})

distance_cys_99_195_r_wt_df = pd.DataFrame({'Distance': distance_cys_r_wt[:, 2],
                                            'Homology Model': np.repeat('ROSETTA', len(distance_cys_r_wt[:, 2]))})

distance_99_195_df = pd.concat(
    [distance_cys_99_195_all_wt_df, distance_cys_99_195_it_wt_df, distance_cys_99_195_m_wt_df,
     distance_cys_99_195_r_wt_df], ignore_index=True)

# RICERCA FRAMES COERENTI CON DISTANZE SMFRET
distance_cys_wt_df_summary = pd.DataFrame({'Distance 10-99': distance_10_99_df['Distance'].loc[distance_10_99_df['Homology Model'] == 'All wt'],
                                           'Distance 10-195': distance_10_195_df['Distance'].loc[distance_10_99_df['Homology Model'] == 'All wt'],
                                           'Distance 99-195': distance_99_195_df['Distance'].loc[distance_10_99_df['Homology Model'] == 'All wt'],
                                           'Homology Model': distance_10_99_df['Homology Model'].loc[distance_10_99_df['Homology Model'] == 'All wt']})

peak_10_99 = 3.78
peak_10_195 = 4.48
peak_99_195 = 4.94

peak_smfret = [peak_10_99, peak_10_195, peak_99_195]

# frames_smfret = distance_cys_wt_df_summary.loc[
#     (distance_cys_wt_df_summary['Distance 10-99'] < 3.88) & (distance_cys_wt_df_summary['Distance 10-99'] > 2.28) &
#     (distance_cys_wt_df_summary['Distance 10-195'] < 4.68) & (distance_cys_wt_df_summary['Distance 10-195'] > 2.98) &
#     (distance_cys_wt_df_summary['Distance 99-195'] < 5.04) & (distance_cys_wt_df_summary['Distance 99-195'] > 3.44)]

dist_frames_peaks = [None for i in range(0, len(distance_cys_wt_df_summary['Distance 10-99']))]

for i in range(0, len(distance_cys_wt_df_summary['Distance 10-99'])):
    dist_frames_peaks[i] = math.dist(peak_smfret, [distance_cys_wt_df_summary.loc[i].at['Distance 10-99'],
                                                   distance_cys_wt_df_summary.loc[i].at['Distance 10-195'],
                                                   distance_cys_wt_df_summary.loc[i].at['Distance 99-195']])

dist_frames_peaks_sorted = sorted(dist_frames_peaks)

mindist = [None for i in range(0, 20)]

for i in range(0, 20):
    mindist[i] = dist_frames_peaks.index(dist_frames_peaks_sorted[i])

mindist_df = distance_cys_wt_df_summary.loc[mindist, ['Distance 10-99', 'Distance 10-195', 'Distance 99-195']]

print(dist_frames_peaks_sorted[0:20])
print(mindist_df)

# SCATTERPLOT 3D

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

p = ax.scatter(distance_cys_wt_df_summary['Distance 10-99'], distance_cys_wt_df_summary['Distance 10-195'],
               distance_cys_wt_df_summary['Distance 99-195'], marker='o', c=dist_frames_peaks, cmap='magma_r')

ax.scatter(peak_10_99, peak_10_195, peak_99_195, marker='o', color='red', s=30)

ax.set_title("Distance from smFRET peaks")
ax.set_xlabel('Distance CYS10-CYS99 (nm)', size=14)
ax.set_ylabel('Distance CYS10-CYS195 (nm)', size=14)
ax.set_zlabel('Distance CYS99-CYS195 (nm)', size=14)

for label in (ax.get_xticklabels() + ax.get_yticklabels()):
    label.set_fontsize(12)

# plt.show()

# fig.savefig('scatter_distance_smfret_wt.png', dpi=320)


# PLOT
sns.set_style('darkgrid')

fig, ax = plt.subplots(1, 1, figsize=(50, 20), gridspec_kw=dict(width_ratios=[4]))

sns.kdeplot(data=distance_10_99_df, x="Distance", hue="Homology Model", linewidth=3, shade=True, common_norm=True)

# FRET PEAKS
plt.axvline(3.78, color='red', linestyle='--', linewidth=3)
# plt.axvline(4.94, color='red', linestyle='--', linewidth=3)

plt.text(3.38, 0.20, 'Peak 2', color='red', fontsize=20)
# plt.text(4.34, 0.20, 'Peak 3', color='red', fontsize=20)

# ENCODERMAP
for i in range(0, len(distance_cys_emap_wt[:, 0])):
    plt.axvline(distance_cys_emap_wt[i, 0], color='blue', linestyle=':', linewidth=2)

# TODO: Calcolare width in nm della gaussiana della E ed inserire la gaussiana nei grafici

sns.set(font_scale=3.5)

ax.set_title("Distance CYS10-CYS99")
ax.set_xlabel('Distance (nm)', size=14)
ax.set_ylabel('Relative frequency', size=14)

for label in (ax.get_xticklabels() + ax.get_yticklabels()):
    label.set_fontsize(12)

# plt.show()

# fig.savefig('distance_10-99_kde_wt.png', dpi=320)
