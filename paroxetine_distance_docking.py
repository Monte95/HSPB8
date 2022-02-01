import mdtraj as md
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import itertools

# DOCKING WT
docking_wt_8 = ["generated_encodermap_wt_8_paroxetine_pose_1.pdb", "generated_encodermap_wt_8_paroxetine_pose_2.pdb",
                "generated_encodermap_wt_8_paroxetine_pose_3.pdb", "generated_encodermap_wt_8_paroxetine_pose_4.pdb",
                "generated_encodermap_wt_8_paroxetine_pose_5.pdb", "generated_encodermap_wt_8_paroxetine_pose_6.pdb",
                "generated_encodermap_wt_8_paroxetine_pose_7.pdb", "generated_encodermap_wt_8_paroxetine_pose_8.pdb",
                "generated_encodermap_wt_8_paroxetine_pose_9.pdb"]

docking_wt_9 = ["generated_encodermap_wt_9_paroxetine_pose_1.pdb", "generated_encodermap_wt_9_paroxetine_pose_2.pdb",
                "generated_encodermap_wt_9_paroxetine_pose_3.pdb", "generated_encodermap_wt_9_paroxetine_pose_4.pdb",
                "generated_encodermap_wt_9_paroxetine_pose_5.pdb", "generated_encodermap_wt_9_paroxetine_pose_6.pdb",
                "generated_encodermap_wt_9_paroxetine_pose_7.pdb", "generated_encodermap_wt_9_paroxetine_pose_8.pdb",
                "generated_encodermap_wt_9_paroxetine_pose_9.pdb"]

docking_wt_10 = ["generated_encodermap_wt_10_paroxetine_pose_1.pdb", "generated_encodermap_wt_10_paroxetine_pose_2.pdb",
                 "generated_encodermap_wt_10_paroxetine_pose_3.pdb", "generated_encodermap_wt_10_paroxetine_pose_4.pdb",
                 "generated_encodermap_wt_10_paroxetine_pose_5.pdb", "generated_encodermap_wt_10_paroxetine_pose_6.pdb",
                 "generated_encodermap_wt_10_paroxetine_pose_7.pdb", "generated_encodermap_wt_10_paroxetine_pose_8.pdb",
                 "generated_encodermap_wt_10_paroxetine_pose_9.pdb"]

docking_wt_14 = ["generated_encodermap_wt_14_paroxetine_pose_1.pdb", "generated_encodermap_wt_14_paroxetine_pose_2.pdb",
                 "generated_encodermap_wt_14_paroxetine_pose_3.pdb", "generated_encodermap_wt_14_paroxetine_pose_4.pdb",
                 "generated_encodermap_wt_14_paroxetine_pose_5.pdb", "generated_encodermap_wt_14_paroxetine_pose_6.pdb",
                 "generated_encodermap_wt_14_paroxetine_pose_7.pdb", "generated_encodermap_wt_14_paroxetine_pose_8.pdb",
                 "generated_encodermap_wt_14_paroxetine_pose_9.pdb"]

# DOCKING K141E
docking_k_5 = ["generated_encodermap_k_5_paroxetine_pose_1.pdb", "generated_encodermap_k_5_paroxetine_pose_2.pdb",
               "generated_encodermap_k_5_paroxetine_pose_3.pdb", "generated_encodermap_k_5_paroxetine_pose_4.pdb",
               "generated_encodermap_k_5_paroxetine_pose_5.pdb", "generated_encodermap_k_5_paroxetine_pose_6.pdb",
               "generated_encodermap_k_5_paroxetine_pose_7.pdb", "generated_encodermap_k_5_paroxetine_pose_8.pdb",
               "generated_encodermap_k_5_paroxetine_pose_9.pdb"]

docking_k_6 = ["generated_encodermap_k_6_paroxetine_pose_1.pdb", "generated_encodermap_k_6_paroxetine_pose_2.pdb",
               "generated_encodermap_k_6_paroxetine_pose_3.pdb", "generated_encodermap_k_6_paroxetine_pose_4.pdb",
               "generated_encodermap_k_6_paroxetine_pose_5.pdb", "generated_encodermap_k_6_paroxetine_pose_6.pdb",
               "generated_encodermap_k_6_paroxetine_pose_7.pdb", "generated_encodermap_k_6_paroxetine_pose_8.pdb",
               "generated_encodermap_k_6_paroxetine_pose_9.pdb"]

docking_k_21 = ["generated_encodermap_k_21_paroxetine_pose_1.pdb", "generated_encodermap_k_21_paroxetine_pose_2.pdb",
                "generated_encodermap_k_21_paroxetine_pose_3.pdb", "generated_encodermap_k_21_paroxetine_pose_4.pdb",
                "generated_encodermap_k_21_paroxetine_pose_5.pdb", "generated_encodermap_k_21_paroxetine_pose_6.pdb",
                "generated_encodermap_k_21_paroxetine_pose_7.pdb", "generated_encodermap_k_21_paroxetine_pose_8.pdb",
                "generated_encodermap_k_21_paroxetine_pose_9.pdb"]

docking_k_22 = ["generated_encodermap_k_22_paroxetine_pose_1.pdb", "generated_encodermap_k_22_paroxetine_pose_2.pdb",
                "generated_encodermap_k_22_paroxetine_pose_3.pdb", "generated_encodermap_k_22_paroxetine_pose_4.pdb",
                "generated_encodermap_k_22_paroxetine_pose_5.pdb", "generated_encodermap_k_22_paroxetine_pose_6.pdb",
                "generated_encodermap_k_22_paroxetine_pose_7.pdb", "generated_encodermap_k_22_paroxetine_pose_8.pdb",
                "generated_encodermap_k_22_paroxetine_pose_9.pdb"]

t_wt_8 = md.load(docking_wt_8)
t_wt_9 = md.load(docking_wt_9)
t_wt_10 = md.load(docking_wt_10)
t_wt_14 = md.load(docking_wt_14)

t_k_5 = md.load(docking_k_5)
t_k_6 = md.load(docking_k_6)
t_k_21 = md.load(docking_k_21)
t_k_22 = md.load(docking_k_22)


paroxetine = [196]
residues = [i for i in range(0, 196)]
residues_string = ["1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18",
                   "19", "20", "21", "22", "23", "24", "25", "26", "27", "28", "29", "30", "31", "32", "33", "34", "35",
                   "36", "37", "38", "39", "40", "41", "42", "43", "44", "45", "46", "47", "48", "49", "50", "51", "52",
                   "53", "54", "55", "56", "57", "58", "59", "60", "61", "62", "63", "64", "65", "66", "67", "68", "69",
                   "70", "71", "72", "73", "74", "75", "76", "77", "78", "79", "80", "81", "82", "83", "84", "85", "86",
                   "87", "88", "89", "90", "91", "92", "93", "94", "95", "96", "97", "98", "99", "100", "101", "102",
                   "103", "104", "105", "106", "107", "108", "109", "110", "111", "112", "113", "114", "115", "116",
                   "117", "118", "119", "120", "121", "122", "123", "124", "125", "126", "127", "128", "129", "130",
                   "131", "132", "133", "134", "135", "136", "137", "138", "139", "140", "141", "142", "143", "144",
                   "145", "146", "147", "148", "149", "150", "151", "152", "153", "154", "155", "156", "157", "158",
                   "159", "160", "161", "162", "163", "164", "165", "166", "167", "168", "169", "170", "171", "172",
                   "173", "174", "175", "176", "177", "178", "179", "180", "181", "182", "183", "184", "185", "186",
                   "187", "188", "189", "190", "191", "192", "193", "194", "195", "196"]

contact_pairs = list(itertools.product(paroxetine, residues))

# COMPUTE CONTACTS
contacts_paroxetine_wt_8, paroxetine_pairs = md.compute_contacts(t_wt_8, contact_pairs)
contacts_paroxetine_wt_9, paroxetine_pairs = md.compute_contacts(t_wt_9, contact_pairs)
contacts_paroxetine_wt_10, paroxetine_pairs = md.compute_contacts(t_wt_10, contact_pairs)
contacts_paroxetine_wt_14, paroxetine_pairs = md.compute_contacts(t_wt_14, contact_pairs)

contacts_paroxetine_k_5, paroxetine_pairs = md.compute_contacts(t_k_5, contact_pairs)
contacts_paroxetine_k_6, paroxetine_pairs = md.compute_contacts(t_k_6, contact_pairs)
contacts_paroxetine_k_21, paroxetine_pairs = md.compute_contacts(t_k_21, contact_pairs)
contacts_paroxetine_k_22, paroxetine_pairs = md.compute_contacts(t_k_22, contact_pairs)

# CONTACTS DATAFRAME
contacts_paroxetine_wt_8_df = pd.DataFrame(contacts_paroxetine_wt_8, columns=residues_string)

contacts_paroxetine_wt_8_df.loc[contacts_paroxetine_wt_8_df["1"] > 0.35] = False
contacts_paroxetine_wt_8_df.loc[contacts_paroxetine_wt_8_df["1"] < 0.35] = True
print(contacts_paroxetine_wt_8_df["1"])

for col in residues_string:
    contacts_paroxetine_wt_8_df.loc[contacts_paroxetine_wt_8_df[col] > 0.35] = False
    contacts_paroxetine_wt_8_df.loc[contacts_paroxetine_wt_8_df[col] < 0.35] = True
print(contacts_paroxetine_wt_8_df)

contacts_paroxetine_wt_8_df_any = contacts_paroxetine_wt_8_df.any()

print(contacts_paroxetine_wt_8_df_any)
contacts_paroxetine_wt_9_df = pd.DataFrame(contacts_paroxetine_wt_9)
for col in contacts_paroxetine_wt_9_df.columns:
    contacts_paroxetine_wt_9_df.loc[contacts_paroxetine_wt_9_df[col] > 3.5, col] = False
    contacts_paroxetine_wt_9_df.loc[contacts_paroxetine_wt_9_df[col] <= 3.5, col] = True

contacts_paroxetine_wt_10_df = pd.DataFrame(contacts_paroxetine_wt_10)
for col in contacts_paroxetine_wt_10_df.columns:
    contacts_paroxetine_wt_10_df.loc[contacts_paroxetine_wt_10_df[col] > 3.5, col] = False
    contacts_paroxetine_wt_10_df.loc[contacts_paroxetine_wt_10_df[col] <= 3.5, col] = True

contacts_paroxetine_wt_14_df = pd.DataFrame(contacts_paroxetine_wt_14)
for col in contacts_paroxetine_wt_14_df.columns:
    contacts_paroxetine_wt_14_df.loc[contacts_paroxetine_wt_14_df[col] > 3.5, col] = False
    contacts_paroxetine_wt_14_df.loc[contacts_paroxetine_wt_14_df[col] <= 3.5, col] = True

contacts_paroxetine_k_5_df = pd.DataFrame(contacts_paroxetine_k_5)
for col in contacts_paroxetine_k_5_df.columns:
    contacts_paroxetine_k_5_df.loc[contacts_paroxetine_k_5_df[col] > 3.5, col] = False
    contacts_paroxetine_k_5_df.loc[contacts_paroxetine_k_5_df[col] <= 3.5, col] = True

contacts_paroxetine_k_6_df = pd.DataFrame(contacts_paroxetine_k_6)
for col in contacts_paroxetine_k_6_df.columns:
    contacts_paroxetine_k_6_df.loc[contacts_paroxetine_k_6_df[col] > 3.5, col] = False
    contacts_paroxetine_k_6_df.loc[contacts_paroxetine_k_6_df[col] <= 3.5, col] = True

contacts_paroxetine_k_21_df = pd.DataFrame(contacts_paroxetine_k_21)
for col in contacts_paroxetine_k_21_df.columns:
    contacts_paroxetine_k_21_df.loc[contacts_paroxetine_k_21_df[col] > 3.5, col] = False
    contacts_paroxetine_k_21_df.loc[contacts_paroxetine_k_21_df[col] <= 3.5, col] = True

contacts_paroxetine_k_22_df = pd.DataFrame(contacts_paroxetine_k_22)
for col in contacts_paroxetine_k_22_df.columns:
    contacts_paroxetine_k_22_df.loc[contacts_paroxetine_k_22_df[col] > 3.5, col] = False
    contacts_paroxetine_k_22_df.loc[contacts_paroxetine_k_22_df[col] <= 3.5, col] = True

# DATAFRAME
paroxetine_distance_wt_8_df = pd.DataFrame({'Paroxetine': [i + 1 for i in paroxetine_pairs[:, 0]],
                                            'Residues': [i + 1 for i in paroxetine_pairs[:, 1]],
                                            'Distance': contacts_paroxetine_wt_8_df_median})

paroxetine_distance_wt_9_df = pd.DataFrame({'Paroxetine': [i + 1 for i in paroxetine_pairs[:, 0]],
                                            'Residues': [i + 1 for i in paroxetine_pairs[:, 1]],
                                            'Distance': contacts_paroxetine_wt_9_df_median})

paroxetine_distance_wt_10_df = pd.DataFrame({'Paroxetine': [i + 1 for i in paroxetine_pairs[:, 0]],
                                             'Residues': [i + 1 for i in paroxetine_pairs[:, 1]],
                                             'Distance': contacts_paroxetine_wt_10_df_median})

paroxetine_distance_wt_14_df = pd.DataFrame({'Paroxetine': [i + 1 for i in paroxetine_pairs[:, 0]],
                                             'Residues': [i + 1 for i in paroxetine_pairs[:, 1]],
                                             'Distance': contacts_paroxetine_wt_14_df_median})

paroxetine_distance_k_5_df = pd.DataFrame({'Paroxetine': [i + 1 for i in paroxetine_pairs[:, 0]],
                                           'Residues': [i + 1 for i in paroxetine_pairs[:, 1]],
                                           'Distance': contacts_paroxetine_k_5_df_median})

paroxetine_distance_k_6_df = pd.DataFrame({'Paroxetine': [i + 1 for i in paroxetine_pairs[:, 0]],
                                           'Residues': [i + 1 for i in paroxetine_pairs[:, 1]],
                                           'Distance': contacts_paroxetine_k_6_df_median})

paroxetine_distance_k_21_df = pd.DataFrame({'Paroxetine': [i + 1 for i in paroxetine_pairs[:, 0]],
                                            'Residues': [i + 1 for i in paroxetine_pairs[:, 1]],
                                            'Distance': contacts_paroxetine_k_21_df_median})

paroxetine_distance_k_22_df = pd.DataFrame({'Paroxetine': [i + 1 for i in paroxetine_pairs[:, 0]],
                                            'Residues': [i + 1 for i in paroxetine_pairs[:, 1]],
                                            'Distance': contacts_paroxetine_k_22_df_median})

# PIVOT
paroxetine_contacts_matrix_wt_8 = paroxetine_distance_wt_8_df.pivot(index='Paroxetine', columns='Residues',
                                                                    values='Distance')

paroxetine_contacts_matrix_wt_9 = paroxetine_distance_wt_9_df.pivot(index='Paroxetine', columns='Residues',
                                                                    values='Distance')

paroxetine_contacts_matrix_wt_10 = paroxetine_distance_wt_10_df.pivot(index='Paroxetine', columns='Residues',
                                                                      values='Distance')

paroxetine_contacts_matrix_wt_14 = paroxetine_distance_wt_14_df.pivot(index='Paroxetine', columns='Residues',
                                                                      values='Distance')

paroxetine_contacts_matrix_k_5 = paroxetine_distance_k_5_df.pivot(index='Paroxetine', columns='Residues',
                                                                  values='Distance')

paroxetine_contacts_matrix_k_6 = paroxetine_distance_k_6_df.pivot(index='Paroxetine', columns='Residues',
                                                                  values='Distance')

paroxetine_contacts_matrix_k_21 = paroxetine_distance_k_21_df.pivot(index='Paroxetine', columns='Residues',
                                                                    values='Distance')

paroxetine_contacts_matrix_k_22 = paroxetine_distance_k_22_df.pivot(index='Paroxetine', columns='Residues',
                                                                    values='Distance')

# vmin = min(min(contacts_paroxetine_wt_df_median), min(contacts_paroxetine_k_df_median))
# vmax = max(max(contacts_paroxetine_wt_df_median), max(contacts_paroxetine_k_df_median))

# HEATMAPS
fig, ax = plt.subplots(4, 2, sharex=True)

sns.heatmap(paroxetine_contacts_matrix_wt_8, cmap="Blues", vmin=0, vmax=5, ax=ax[0, 0], cbar_kws={'label': 'Distance (nm)'})
sns.heatmap(paroxetine_contacts_matrix_wt_9, cmap="Blues", vmin=0, vmax=5, ax=ax[1, 0], cbar_kws={'label': 'Distance (nm)'})
sns.heatmap(paroxetine_contacts_matrix_wt_10, cmap="Blues", vmin=0, vmax=5, ax=ax[2, 0], cbar_kws={'label': 'Distance (nm)'})
sns.heatmap(paroxetine_contacts_matrix_wt_14, cmap="Blues", vmin=0, vmax=5, ax=ax[3, 0], cbar_kws={'label': 'Distance (nm)'})

sns.heatmap(paroxetine_contacts_matrix_k_5, cmap="Blues", vmin=0, vmax=5, ax=ax[0, 1], cbar_kws={'label': 'Distance (nm)'})
sns.heatmap(paroxetine_contacts_matrix_k_6, cmap="Blues", vmin=0, vmax=5, ax=ax[1, 1], cbar_kws={'label': 'Distance (nm)'})
sns.heatmap(paroxetine_contacts_matrix_k_21, cmap="Blues", vmin=0, vmax=5, ax=ax[2, 1], cbar_kws={'label': 'Distance (nm)'})
sns.heatmap(paroxetine_contacts_matrix_k_22, cmap="Blues", vmin=0, vmax=5, ax=ax[3, 1], cbar_kws={'label': 'Distance (nm)'})

# ax.set_title("wt", fontsize=20)
# ax[1].set_title("K141E", fontsize=20)

# ax.set_xlabel(None)
# ax[1].set_xlabel('Residue', fontsize=20)

# ax.set_ylabel('Paroxetine', fontsize=20)
# ax[1].set_ylabel('Paroxetine', fontsize=20)

plt.show()

# fig.savefig("paroxetine_contactmaps_median_wt_k.png", dpi=320)








