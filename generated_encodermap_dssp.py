import mdtraj as md
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

NUM_RESIDUES = 196

# t = md.load("generated_encodermap_wt_19.pdb")

# dssp = np.transpose(md.compute_dssp(t, simplified=False))
# dssp_df = pd.DataFrame(dssp)

# secstr_tot_dssp = (dssp_df.apply(pd.value_counts) / NUM_RESIDUES).fillna(0)


percentage = np.array([0.061, 0.036, 0.036, 0.041, 0.102, 0.122, 0.107, 0.092, 0.061, 0.026, 0.061, 0.071, 0.087, 0.046,
                       0, 0.041, 0.056, 0.077, 0.077, 0.051, 0.056, 0.077, 0.061, 0.056, # Helix
                       0.031, 0.077, 0.046, 0.046, 0, 0.015, 0.015, 0.015, 0.061, 0.061, 0.020, 0.015, 0.015, 0.015,
                       0.015, 0, 0, 0.046, 0.061, 0.061, 0.061, 0.061, 0.020, 0.0446, # Distorted Helix
                       0.179, 0.179, 0.230, 0.179, 0.158, 0.158, 0.199, 0.204, 0.117, 0.173, 0.204, 0.199, 0.194,
                       0.224, 0.245, 0.245, 0.163, 0.168, 0.168, 0.168, 0.132, 0.153, 0.179, 0.168, # Strand
                       0.010, 0.020, 0, 0.020, 0.031, 0.020, 0.020, 0.010, 0.030, 0.010, 0.020, 0.010, 0.026, 0, 0.020,
                       0.005, 0.031, 0.010, 0, 0.020, 0.010, 0.031, 0.020, 0.010, # Distorted Strand
                       0.122, 0.097, 0.153, 0.194, 0.091, 0.081, 0.112, 0.133, 0.153, 0.163, 0.107, 0.112, 0.102, 0.209,
                       0.143, 0.153, 0.163, 0.112, 0.091, 0.138, 0.128, 0.128, 0.112, 0.138, # Turn
                       0.597, 0.592, 0.536, 0.520, 0.617, 0.602, 0.546, 0.546, 0.577, 0.567, 0.587, 0.592, 0.576, 0.505,
                       0.576, 0.556, 0.587, 0.586, 0.602, 0.562, 0.612, 0.551, 0.577, 0.576 # Unordered
                       ]) * 100

secondary_structure = ["Helix", "Helix", "Helix", "Helix", "Helix", "Helix", "Helix", "Helix", "Helix", "Helix",
                       "Helix", "Helix", "Helix", "Helix", "Helix", "Helix", "Helix", "Helix", "Helix", "Helix",
                       "Helix", "Helix", "Helix", "Helix", #
                       "Distorted Helix", "Distorted Helix", "Distorted Helix", "Distorted Helix", "Distorted Helix",
                       "Distorted Helix", "Distorted Helix", "Distorted Helix", "Distorted Helix", "Distorted Helix",
                       "Distorted Helix", "Distorted Helix", "Distorted Helix", "Distorted Helix", "Distorted Helix",
                       "Distorted Helix", "Distorted Helix", "Distorted Helix", "Distorted Helix", "Distorted Helix",
                       "Distorted Helix", "Distorted Helix", "Distorted Helix", "Distorted Helix", #
                       "Strand", "Strand", "Strand", "Strand", "Strand", "Strand", "Strand", "Strand", "Strand",
                       "Strand", "Strand", "Strand", "Strand", "Strand", "Strand", "Strand", "Strand", "Strand",
                       "Strand", "Strand", "Strand", "Strand", "Strand", "Strand", #
                       "Distorted Strand", "Distorted Strand", "Distorted Strand", "Distorted Strand",
                       "Distorted Strand", "Distorted Strand", "Distorted Strand", "Distorted Strand",
                       "Distorted Strand", "Distorted Strand", "Distorted Strand", "Distorted Strand",
                       "Distorted Strand", "Distorted Strand", "Distorted Strand", "Distorted Strand",
                       "Distorted Strand", "Distorted Strand", "Distorted Strand", "Distorted Strand",
                       "Distorted Strand", "Distorted Strand", "Distorted Strand", "Distorted Strand", #
                       "Turn", "Turn", "Turn", "Turn", "Turn", "Turn", "Turn", "Turn", "Turn", "Turn", "Turn", "Turn",
                       "Turn", "Turn", "Turn", "Turn", "Turn", "Turn", "Turn", "Turn", "Turn", "Turn", "Turn", "Turn",#
                       "Unordered", "Unordered", "Unordered", "Unordered", "Unordered", "Unordered", "Unordered",
                       "Unordered", "Unordered", "Unordered", "Unordered", "Unordered", "Unordered", "Unordered",
                       "Unordered", "Unordered", "Unordered", "Unordered", "Unordered", "Unordered", "Unordered",
                       "Unordered", "Unordered", "Unordered"]

variant = ["wt 1", "wt 2", "wt 3", "wt 4", "wt 5", "wt 6", "wt 7", "wt 8", "wt 9", "wt 10", "wt 11", "wt 12", "wt 13",
           "K141E 1", "K141E 2", "K141E 3", "K141E 4", "K141E 5", "K141E 6", "K141E 7", "K141E 8", "K141E 9",
           "Median wt", "Median K141E", #
           "wt 1", "wt 2", "wt 3", "wt 4", "wt 5", "wt 6", "wt 7", "wt 8", "wt 9", "wt 10", "wt 11", "wt 12", "wt 13",
           "K141E 1", "K141E 2", "K141E 3", "K141E 4", "K141E 5", "K141E 6", "K141E 7", "K141E 8", "K141E 9",
           "Median wt", "Median K141E", #
           "wt 1", "wt 2", "wt 3", "wt 4", "wt 5", "wt 6", "wt 7", "wt 8", "wt 9", "wt 10", "wt 11", "wt 12", "wt 13",
           "K141E 1", "K141E 2", "K141E 3", "K141E 4", "K141E 5", "K141E 6", "K141E 7", "K141E 8", "K141E 9",
           "Median wt", "Median K141E", #
           "wt 1", "wt 2", "wt 3", "wt 4", "wt 5", "wt 6", "wt 7", "wt 8", "wt 9", "wt 10", "wt 11", "wt 12", "wt 13",
           "K141E 1", "K141E 2", "K141E 3", "K141E 4", "K141E 5", "K141E 6", "K141E 7", "K141E 8", "K141E 9",
           "Median wt", "Median K141E", #
           "wt 1", "wt 2", "wt 3", "wt 4", "wt 5", "wt 6", "wt 7", "wt 8", "wt 9", "wt 10", "wt 11", "wt 12", "wt 13",
           "K141E 1", "K141E 2", "K141E 3", "K141E 4", "K141E 5", "K141E 6", "K141E 7", "K141E 8", "K141E 9",
           "Median wt", "Median K141E", #
           "wt 1", "wt 2", "wt 3", "wt 4", "wt 5", "wt 6", "wt 7", "wt 8", "wt 9", "wt 10", "wt 11", "wt 12", "wt 13",
           "K141E 1", "K141E 2", "K141E 3", "K141E 4", "K141E 5", "K141E 6", "K141E 7", "K141E 8", "K141E 9",
           "Median wt", "Median K141E"]

source = ["Encodermap", "Encodermap", "Encodermap", "Encodermap", "Encodermap", "Encodermap", "Encodermap",
          "Encodermap", "Encodermap", "Encodermap", "Encodermap", "Encodermap", "Encodermap", "Encodermap",
          "Encodermap", "Encodermap", "Encodermap", "Encodermap", "Encodermap", "Encodermap", "Encodermap",
          "Encodermap", "Encodermap", "Encodermap", #
          "Encodermap", "Encodermap", "Encodermap", "Encodermap", "Encodermap", "Encodermap", "Encodermap",
          "Encodermap", "Encodermap", "Encodermap", "Encodermap", "Encodermap", "Encodermap", "Encodermap",
          "Encodermap", "Encodermap", "Encodermap", "Encodermap", "Encodermap", "Encodermap", "Encodermap",
          "Encodermap", "Encodermap", "Encodermap", #
          "Encodermap", "Encodermap", "Encodermap", "Encodermap", "Encodermap", "Encodermap", "Encodermap",
          "Encodermap", "Encodermap", "Encodermap", "Encodermap", "Encodermap", "Encodermap", "Encodermap",
          "Encodermap", "Encodermap", "Encodermap", "Encodermap", "Encodermap", "Encodermap", "Encodermap",
          "Encodermap", "Encodermap", "Encodermap", #
          "Encodermap", "Encodermap", "Encodermap", "Encodermap", "Encodermap", "Encodermap", "Encodermap",
          "Encodermap", "Encodermap", "Encodermap", "Encodermap", "Encodermap", "Encodermap", "Encodermap",
          "Encodermap", "Encodermap", "Encodermap", "Encodermap", "Encodermap", "Encodermap", "Encodermap",
          "Encodermap", "Encodermap", "Encodermap", #
          "Encodermap", "Encodermap", "Encodermap", "Encodermap", "Encodermap", "Encodermap", "Encodermap",
          "Encodermap", "Encodermap", "Encodermap", "Encodermap", "Encodermap", "Encodermap", "Encodermap",
          "Encodermap", "Encodermap", "Encodermap", "Encodermap", "Encodermap", "Encodermap", "Encodermap",
          "Encodermap", "Encodermap", "Encodermap", #
          "Encodermap", "Encodermap", "Encodermap", "Encodermap", "Encodermap", "Encodermap", "Encodermap",
          "Encodermap", "Encodermap", "Encodermap", "Encodermap", "Encodermap", "Encodermap", "Encodermap",
          "Encodermap", "Encodermap", "Encodermap", "Encodermap", "Encodermap", "Encodermap", "Encodermap",
          "Encodermap", "Encodermap", "Encodermap"]


dssp_distorted_df = pd.DataFrame(
    {"Percentage": percentage,
     "Secondary structure": secondary_structure,
     "Variant": variant,
     "Source": source})

########
# PLOT
########

colors = ['Blue', 'Blue', 'Blue', 'Blue', 'Blue', 'Blue', 'Blue', 'Blue', 'Blue', 'Blue', 'Blue', 'Blue', 'Blue',
          'Orange', 'Orange', 'Orange', 'Orange', 'Orange', 'Orange', 'Orange', 'Orange', 'Orange', "Blue", "Orange"]

sns.set_style('darkgrid')

g = sns.catplot(data=dssp_distorted_df, x='Secondary structure', y='Percentage', hue='Variant',
                kind='bar', legend_out=False, palette=colors, legend=False)

ax = g.facet_axis(0, 0)

i = 0
for c in ax.containers:
    if i in [22, 23]:
        labels = [f'{(v.get_height()):.1f}' for v in c]
        ax.bar_label(c, labels=labels, label_type='edge', fontsize=14)


for label in (ax.get_xticklabels() + ax.get_yticklabels()):
    label.set_fontsize(14)


ax.set_xlabel('Secondary structure', size=14)
ax.set_ylabel('Percentage', size=14)

ax.set_title("Encodermap", size=14)


plt.show()

#plt.savefig("secstr_comparison_CD_dssp.png", dpi=320)





