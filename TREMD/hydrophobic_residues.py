import mdtraj as md
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import itertools

ref_structure_wt = "it_wt.gro"
trajectory_wt = "all_wt.xtc"

ref_structure_k = "it_k.gro"
trajectory_k = "all_k.xtc"

t_wt = md.load(trajectory_wt, top=ref_structure_wt)
t_k = md.load(trajectory_k, top=ref_structure_k)


def read_sasa_data(document_name):
  x, y = [], []
  with open(document_name) as f:
    for line in f:
        cols = line.split()

        x.append(float(cols[0]))
        y.append(float(cols[1]))
  return x, y


hydrophobic_resides_sequence = [1, 2, 6, 7, 8, 13, 16, 20, 21, 25, 26, 30, 31, 35, 37,
                                39, 40, 41, 44, 46, 48, 49, 51, 52, 53, 54, 56, 59, 60, 61, 64, 68,
                                69, 70, 73, 75, 77, 79, 81, 82, 83, 88, 89, 90, 91, 92,
                                95, 96, 98, 100, 102, 105, 107, 110, 111, 112, 119, 121, 134,
                                135, 139, 143, 145, 146, 147, 149, 151, 152, 154, 155, 156,
                                158, 160, 163, 164, 165, 166, 168, 169, 171, 172, 173, 177,
                                182, 186, 187, 193]

hydrophobic_resides_id = [i - 1 for i in hydrophobic_resides_sequence]

hydrophobic_resides_acd = [95, 97, 99, 101, 104, 106, 109, 110, 111, 118, 120]

# contact_pairs = list(itertools.product(hydrophobic_resides_acd, hydrophobic_resides_id))
# contacts_hydrophobic_wt, hydrophobic_pairs = md.compute_contacts(t_wt, contact_pairs, scheme='sidechain')
#
# contact_pairs = list(itertools.product(hydrophobic_resides_acd, hydrophobic_resides_id))
# contacts_hydrophobic_k, hydrophobic_pairs = md.compute_contacts(t_k, contact_pairs, scheme='sidechain')
#
# contacts_hydrophobic_wt_df = pd.DataFrame(contacts_hydrophobic_wt)
# contacts_hydrophobic_wt_df_median = contacts_hydrophobic_wt_df.median(axis=0)
#
# contacts_hydrophobic_k_df = pd.DataFrame(contacts_hydrophobic_k)
# contacts_hydrophobic_k_df_median = contacts_hydrophobic_k_df.median(axis=0)
#
# hydrophobic_distance_wt_df = pd.DataFrame({'Hydrophobic residue': [i + 1 for i in hydrophobic_pairs[:, 0]],
#                                            'Residue': [i + 1 for i in hydrophobic_pairs[:, 1]],
#                                            'Distance': contacts_hydrophobic_wt_df_median})
#
# hydrophobic_distance_k_df = pd.DataFrame({'Hydrophobic residue': [i + 1 for i in hydrophobic_pairs[:, 0]],
#                                            'Residue': [i + 1 for i in hydrophobic_pairs[:, 1]],
#                                            'Distance': contacts_hydrophobic_k_df_median})
#
# hydrophobic_contacts_matrix_wt = hydrophobic_distance_wt_df.pivot(index='Hydrophobic residue', columns='Residue',
#                                                                   values='Distance')
#
# hydrophobic_contacts_matrix_k = hydrophobic_distance_k_df.pivot(index='Hydrophobic residue', columns='Residue',
#                                                                   values='Distance')
#
# vmin = min(min(contacts_hydrophobic_wt_df_median), min(contacts_hydrophobic_k_df_median))
# vmax = max(max(contacts_hydrophobic_wt_df_median), max(contacts_hydrophobic_k_df_median))
#
# fig, ax = plt.subplots(2, 1, sharex=True)
#
# sns.heatmap(hydrophobic_contacts_matrix_wt, cmap="rocket_r", vmin=vmin, vmax=vmax, ax=ax[0], cbar_kws={'label': 'Distance (nm)'})
# sns.heatmap(hydrophobic_contacts_matrix_k, cmap="rocket_r", vmin=vmin, vmax=vmax, ax=ax[1], cbar_kws={'label': 'Distance (nm)'})
#
# ax[0].set_title("wt", fontsize=20)
# ax[1].set_title("K141E", fontsize=20)
#
# ax[0].set_xlabel(None)
# ax[1].set_xlabel('Residue', fontsize=20)
#
# ax[0].set_ylabel('ACD hydrophobic residue', fontsize=20)
# ax[1].set_ylabel('ACD hydrophobic residue', fontsize=20)
#
# plt.show()
#
# fig.savefig("hydrophobic_ACD_contactmaps_median_wt_k.png", dpi=320)

# SASA

t, sasa_acd_96_wt = read_sasa_data("sasa_tot_acd_wt_96.xvg")
t, sasa_acd_98_wt = read_sasa_data("sasa_tot_acd_wt_98.xvg")
t, sasa_acd_100_wt = read_sasa_data("sasa_tot_acd_wt_100.xvg")
t, sasa_acd_102_wt = read_sasa_data("sasa_tot_acd_wt_102.xvg")
t, sasa_acd_105_wt = read_sasa_data("sasa_tot_acd_wt_105.xvg")
t, sasa_acd_107_wt = read_sasa_data("sasa_tot_acd_wt_107.xvg")
t, sasa_acd_110_wt = read_sasa_data("sasa_tot_acd_wt_110.xvg")
t, sasa_acd_111_wt = read_sasa_data("sasa_tot_acd_wt_111.xvg")
t, sasa_acd_112_wt = read_sasa_data("sasa_tot_acd_wt_112.xvg")
t, sasa_acd_119_wt = read_sasa_data("sasa_tot_acd_wt_119.xvg")
t, sasa_acd_121_wt = read_sasa_data("sasa_tot_acd_wt_121.xvg")

t, sasa_acd_96_k = read_sasa_data("sasa_tot_acd_k_96.xvg")
t, sasa_acd_98_k = read_sasa_data("sasa_tot_acd_k_98.xvg")
t, sasa_acd_100_k = read_sasa_data("sasa_tot_acd_k_100.xvg")
t, sasa_acd_102_k = read_sasa_data("sasa_tot_acd_k_102.xvg")
t, sasa_acd_105_k = read_sasa_data("sasa_tot_acd_k_105.xvg")
t, sasa_acd_107_k = read_sasa_data("sasa_tot_acd_k_107.xvg")
t, sasa_acd_110_k = read_sasa_data("sasa_tot_acd_k_110.xvg")
t, sasa_acd_111_k = read_sasa_data("sasa_tot_acd_k_111.xvg")
t, sasa_acd_112_k = read_sasa_data("sasa_tot_acd_k_112.xvg")
t, sasa_acd_119_k = read_sasa_data("sasa_tot_acd_k_119.xvg")
t, sasa_acd_121_k = read_sasa_data("sasa_tot_acd_k_121.xvg")

t, sasa_tot_wt = read_sasa_data("sasa_tot_all_wt.xvg")

t, sasa_tot_k = read_sasa_data("sasa_tot_all_k.xvg")

# DATAFRAME WT
acd_96_wt_sasa_df = pd.DataFrame({'Time': t,
                                  'SASA': sasa_tot_wt,
                                  'SASA Residue': sasa_acd_96_wt,
                                  'SASA Residue/SASA': np.array(sasa_acd_96_wt) / np.array(sasa_tot_wt),
                                  'Variant': np.repeat('wt', len(sasa_tot_wt)),
                                  'Residue': np.repeat('96', len(sasa_tot_wt))})

acd_98_wt_sasa_df = pd.DataFrame({'Time': t,
                                  'SASA': sasa_tot_wt,
                                  'SASA Residue': sasa_acd_98_wt,
                                  'SASA Residue/SASA': np.array(sasa_acd_98_wt) / np.array(sasa_tot_wt),
                                  'Variant': np.repeat('wt', len(sasa_tot_wt)),
                                  'Residue': np.repeat('98', len(sasa_tot_wt))})

acd_100_wt_sasa_df = pd.DataFrame({'Time': t,
                                  'SASA': sasa_tot_wt,
                                  'SASA Residue': sasa_acd_100_wt,
                                  'SASA Residue/SASA': np.array(sasa_acd_100_wt) / np.array(sasa_tot_wt),
                                  'Variant': np.repeat('wt', len(sasa_tot_wt)),
                                  'Residue': np.repeat('100', len(sasa_tot_wt))})

acd_102_wt_sasa_df = pd.DataFrame({'Time': t,
                                  'SASA': sasa_tot_wt,
                                  'SASA Residue': sasa_acd_102_wt,
                                  'SASA Residue/SASA': np.array(sasa_acd_102_wt) / np.array(sasa_tot_wt),
                                  'Variant': np.repeat('wt', len(sasa_tot_wt)),
                                  'Residue': np.repeat('102', len(sasa_tot_wt))})

acd_105_wt_sasa_df = pd.DataFrame({'Time': t,
                                  'SASA': sasa_tot_wt,
                                  'SASA Residue': sasa_acd_105_wt,
                                  'SASA Residue/SASA': np.array(sasa_acd_105_wt) / np.array(sasa_tot_wt),
                                  'Variant': np.repeat('wt', len(sasa_tot_wt)),
                                  'Residue': np.repeat('105', len(sasa_tot_wt))})

acd_107_wt_sasa_df = pd.DataFrame({'Time': t,
                                  'SASA': sasa_tot_wt,
                                  'SASA Residue': sasa_acd_107_wt,
                                  'SASA Residue/SASA': np.array(sasa_acd_107_wt) / np.array(sasa_tot_wt),
                                  'Variant': np.repeat('wt', len(sasa_tot_wt)),
                                  'Residue': np.repeat('107', len(sasa_tot_wt))})

acd_110_wt_sasa_df = pd.DataFrame({'Time': t,
                                  'SASA': sasa_tot_wt,
                                  'SASA Residue': sasa_acd_110_wt,
                                  'SASA Residue/SASA': np.array(sasa_acd_110_wt) / np.array(sasa_tot_wt),
                                  'Variant': np.repeat('wt', len(sasa_tot_wt)),
                                  'Residue': np.repeat('110', len(sasa_tot_wt))})

acd_111_wt_sasa_df = pd.DataFrame({'Time': t,
                                  'SASA': sasa_tot_wt,
                                  'SASA Residue': sasa_acd_111_wt,
                                  'SASA Residue/SASA': np.array(sasa_acd_111_wt) / np.array(sasa_tot_wt),
                                  'Variant': np.repeat('wt', len(sasa_tot_wt)),
                                  'Residue': np.repeat('111', len(sasa_tot_wt))})

acd_112_wt_sasa_df = pd.DataFrame({'Time': t,
                                  'SASA': sasa_tot_wt,
                                  'SASA Residue': sasa_acd_112_wt,
                                  'SASA Residue/SASA': np.array(sasa_acd_112_wt) / np.array(sasa_tot_wt),
                                  'Variant': np.repeat('wt', len(sasa_tot_wt)),
                                  'Residue': np.repeat('112', len(sasa_tot_wt))})

acd_119_wt_sasa_df = pd.DataFrame({'Time': t,
                                  'SASA': sasa_tot_wt,
                                  'SASA Residue': sasa_acd_119_wt,
                                  'SASA Residue/SASA': np.array(sasa_acd_119_wt) / np.array(sasa_tot_wt),
                                  'Variant': np.repeat('wt', len(sasa_tot_wt)),
                                  'Residue': np.repeat('119', len(sasa_tot_wt))})

acd_121_wt_sasa_df = pd.DataFrame({'Time': t,
                                  'SASA': sasa_tot_wt,
                                  'SASA Residue': sasa_acd_121_wt,
                                  'SASA Residue/SASA': np.array(sasa_acd_121_wt) / np.array(sasa_tot_wt),
                                  'Variant': np.repeat('wt', len(sasa_tot_wt)),
                                  'Residue': np.repeat('121', len(sasa_tot_wt))})

# DATAFRAME K141E

acd_96_k_sasa_df = pd.DataFrame({'Time': t,
                                  'SASA': sasa_tot_k,
                                  'SASA Residue': sasa_acd_96_k,
                                  'SASA Residue/SASA': np.array(sasa_acd_96_k) / np.array(sasa_tot_k),
                                  'Variant': np.repeat('K141E', len(sasa_tot_k)),
                                  'Residue': np.repeat('96', len(sasa_tot_k))})

acd_98_k_sasa_df = pd.DataFrame({'Time': t,
                                  'SASA': sasa_tot_k,
                                  'SASA Residue': sasa_acd_98_k,
                                  'SASA Residue/SASA': np.array(sasa_acd_98_k) / np.array(sasa_tot_k),
                                  'Variant': np.repeat('K141E', len(sasa_tot_k)),
                                  'Residue': np.repeat('98', len(sasa_tot_k))})

acd_100_k_sasa_df = pd.DataFrame({'Time': t,
                                  'SASA': sasa_tot_k,
                                  'SASA Residue': sasa_acd_100_k,
                                  'SASA Residue/SASA': np.array(sasa_acd_100_k) / np.array(sasa_tot_k),
                                  'Variant': np.repeat('K141E', len(sasa_tot_k)),
                                  'Residue': np.repeat('100', len(sasa_tot_k))})

acd_102_k_sasa_df = pd.DataFrame({'Time': t,
                                  'SASA': sasa_tot_k,
                                  'SASA Residue': sasa_acd_102_k,
                                  'SASA Residue/SASA': np.array(sasa_acd_102_k) / np.array(sasa_tot_k),
                                  'Variant': np.repeat('K141E', len(sasa_tot_k)),
                                  'Residue': np.repeat('102', len(sasa_tot_k))})

acd_105_k_sasa_df = pd.DataFrame({'Time': t,
                                  'SASA': sasa_tot_k,
                                  'SASA Residue': sasa_acd_105_k,
                                  'SASA Residue/SASA': np.array(sasa_acd_105_k) / np.array(sasa_tot_k),
                                  'Variant': np.repeat('K141E', len(sasa_tot_k)),
                                  'Residue': np.repeat('105', len(sasa_tot_k))})

acd_107_k_sasa_df = pd.DataFrame({'Time': t,
                                  'SASA': sasa_tot_k,
                                  'SASA Residue': sasa_acd_107_k,
                                  'SASA Residue/SASA': np.array(sasa_acd_107_k) / np.array(sasa_tot_k),
                                  'Variant': np.repeat('K141E', len(sasa_tot_k)),
                                  'Residue': np.repeat('107', len(sasa_tot_k))})

acd_110_k_sasa_df = pd.DataFrame({'Time': t,
                                  'SASA': sasa_tot_k,
                                  'SASA Residue': sasa_acd_110_k,
                                  'SASA Residue/SASA': np.array(sasa_acd_110_k) / np.array(sasa_tot_k),
                                  'Variant': np.repeat('K141E', len(sasa_tot_k)),
                                  'Residue': np.repeat('110', len(sasa_tot_k))})

acd_111_k_sasa_df = pd.DataFrame({'Time': t,
                                  'SASA': sasa_tot_k,
                                  'SASA Residue': sasa_acd_111_k,
                                  'SASA Residue/SASA': np.array(sasa_acd_111_k) / np.array(sasa_tot_k),
                                  'Variant': np.repeat('K141E', len(sasa_tot_k)),
                                  'Residue': np.repeat('111', len(sasa_tot_k))})

acd_112_k_sasa_df = pd.DataFrame({'Time': t,
                                  'SASA': sasa_tot_k,
                                  'SASA Residue': sasa_acd_112_k,
                                  'SASA Residue/SASA': np.array(sasa_acd_112_k) / np.array(sasa_tot_k),
                                  'Variant': np.repeat('K141E', len(sasa_tot_k)),
                                  'Residue': np.repeat('112', len(sasa_tot_k))})

acd_119_k_sasa_df = pd.DataFrame({'Time': t,
                                  'SASA': sasa_tot_k,
                                  'SASA Residue': sasa_acd_119_k,
                                  'SASA Residue/SASA': np.array(sasa_acd_119_k) / np.array(sasa_tot_k),
                                  'Variant': np.repeat('K141E', len(sasa_tot_k)),
                                  'Residue': np.repeat('119', len(sasa_tot_k))})

acd_121_k_sasa_df = pd.DataFrame({'Time': t,
                                  'SASA': sasa_tot_k,
                                  'SASA Residue': sasa_acd_121_k,
                                  'SASA Residue/SASA': np.array(sasa_acd_121_k) / np.array(sasa_tot_k),
                                  'Variant': np.repeat('K141E', len(sasa_tot_k)),
                                  'Residue': np.repeat('121', len(sasa_tot_k))})

acd_singoli_wt_k_sasa_df = pd.concat([acd_96_wt_sasa_df, acd_98_wt_sasa_df, acd_100_wt_sasa_df, acd_102_wt_sasa_df,
                                      acd_105_wt_sasa_df, acd_107_wt_sasa_df, acd_110_wt_sasa_df, acd_111_wt_sasa_df,
                                      acd_112_wt_sasa_df, acd_119_wt_sasa_df, acd_121_wt_sasa_df, acd_96_k_sasa_df,
                                      acd_98_k_sasa_df, acd_100_k_sasa_df, acd_102_k_sasa_df,
                                      acd_105_k_sasa_df, acd_107_k_sasa_df, acd_110_k_sasa_df, acd_111_k_sasa_df,
                                      acd_112_k_sasa_df, acd_119_k_sasa_df, acd_121_k_sasa_df],
                                     ignore_index=True)

# PLOT
sns.set_style('darkgrid')

fig, ax = plt.subplots(1, figsize=(50, 20), gridspec_kw=dict(width_ratios=[4]))

# sns.violinplot(data=acd_singoli_wt_k_sasa_df, x="Residue", y='SASA Residue/SASA', hue="Variant", split=True)

sns.boxplot(data=acd_singoli_wt_k_sasa_df, x='Residue', y='SASA Residue', hue='Variant', linewidth=2.5)

ax.set_title("ACD hSASA", fontsize=20)
ax.set_xlabel('Residue', fontsize=20)
ax.set_ylabel('hSASA', fontsize=20)

for label in (ax.get_xticklabels() + ax.get_yticklabels()):
    label.set_fontsize(14)

leg = ax.legend(prop={"size": 20})


plt.show()

# fig.savefig('boxplot_acd_hsasa_wt_k.png', dpi=320)
