import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

# SENZA DISTORTED
cd_df = pd.DataFrame({"Fraction": [0.052, 0.057, 0.365, 0.320, 0.208, 0.165, 0.374, 0.458],
                      "Secondary structure": ["Helix", "Helix", "Strand", "Strand", "Turn", "Turn", "Unordered",
                                              "Unordered"],
                      "Variant": ["wt", "K141E", "wt", "K141E", "wt", "K141E", "wt", "K141E"],
                      "Source": ["CD", "CD", "CD", "CD", "CD", "CD", "CD", "CD"]})

dssp_df = pd.DataFrame({"Fraction": [0.098, 0.082, 0.202, 0.197, 0.131, 0.141, 0.568, 0.580],
                        "Secondary structure": ["Helix", "Helix", "Strand", "Strand", "Turn", "Turn", "Unordered",
                                                "Unordered"],
                        "Variant": ["wt", "K141E", "wt", "K141E", "wt", "K141E", "wt", "K141E"],
                        "Source": ["Simulations", "Simulations", "Simulations", "Simulations", "Simulations",
                                   "Simulations", "Simulations", "Simulations"]})

secstr_df = pd.concat([cd_df, dssp_df], ignore_index=True)

# CON DISTORTED

cd_distorted_df = pd.DataFrame(
    {"Percentage": np.array([0.012, 0.012, 0.040, 0.045, 0.251, 0.228, 0.114, 0.092, 0.208, 0.165, 0.374, 0.458])*100,
     "Secondary structure": ["Helix", "Helix", "Distorted Helix", "Distorted Helix", "Strand", "Strand",
                             "Distorted Strand", "Distorted Strand", "Turn", "Turn", "Unordered", "Unordered"],
     "Variant": ["wt", "K141E", "wt", "K141E", "wt", "K141E", "wt", "K141E", "wt", "K141E", "wt", "K141E"],
     "Source": ["CD", "CD", "CD", "CD", "CD", "CD", "CD", "CD", "CD", "CD", "CD", "CD"]})

dssp_distorted_df = pd.DataFrame(
    {"Percentage": np.array([0.067, 0.048, 0.031, 0.033, 0.184, 0.183, 0.019, 0.014, 0.131, 0.141, 0.568, 0.580])*100,
     "Secondary structure": ["Helix", "Helix", "Distorted Helix", "Distorted Helix", "Strand", "Strand",
                             "Distorted Strand", "Distorted Strand", "Turn", "Turn", "Unordered", "Unordered"],
     "Variant": ["wt", "K141E", "wt", "K141E", "wt", "K141E", "wt", "K141E", "wt", "K141E", "wt", "K141E"],
     "Source": ["Simulations", "Simulations", "Simulations", "Simulations", "Simulations", "Simulations", "Simulations",
                "Simulations", "Simulations", "Simulations", "Simulations", "Simulations"]})

encodermap_distorted_df = pd.DataFrame(
    {"Percentage": np.array([0.061, 0.056, 0.020, 0.0446, 0.179, 0.168, 0.020, 0.010, 0.112, 0.138, 0.577, 0.576])*100,
     "Secondary structure": ["Helix", "Helix", "Distorted Helix", "Distorted Helix", "Strand", "Strand",
                             "Distorted Strand", "Distorted Strand", "Turn", "Turn", "Unordered", "Unordered"],
     "Variant": ["wt", "K141E", "wt", "K141E", "wt", "K141E", "wt", "K141E", "wt", "K141E", "wt", "K141E"],
     "Source": ["Encodermap", "Encodermap", "Encodermap", "Encodermap", "Encodermap", "Encodermap", "Encodermap",
                "Encodermap", "Encodermap", "Encodermap", "Encodermap", "Encodermap"]})

dssp_distorted_df_no_parox = pd.DataFrame(
    {"Percentage": np.array([0.067, 0.048, 0.031, 0.033, 0.184, 0.183, 0.019, 0.014, 0.131, 0.141, 0.568, 0.580])*100,
     "Secondary structure": ["Helix", "Helix", "Distorted Helix", "Distorted Helix", "Strand", "Strand",
                             "Distorted Strand", "Distorted Strand", "Turn", "Turn", "Unordered", "Unordered"],
     "Variant": ["wt", "K141E", "wt", "K141E", "wt", "K141E", "wt", "K141E", "wt", "K141E", "wt", "K141E"],
     "Source": ["No paroxetine", "No paroxetine", "No paroxetine", "No paroxetine", "No paroxetine", "No paroxetine",
                "No paroxetine", "No paroxetine", "No paroxetine", "No paroxetine", "No paroxetine", "No paroxetine"]})

paroxetine_distorted_df = pd.DataFrame(
    {"Percentage": np.array([0.060, 0.054, 0.039, 0.036, 0.181, 0.189, 0.019, 0.016, 0.137, 0.137, 0.563, 0.566])*100,
     "Secondary structure": ["Helix", "Helix", "Distorted Helix", "Distorted Helix", "Strand", "Strand",
                             "Distorted Strand", "Distorted Strand", "Turn", "Turn", "Unordered", "Unordered"],
     "Variant": ["wt", "K141E", "wt", "K141E", "wt", "K141E", "wt", "K141E", "wt", "K141E", "wt", "K141E"],
     "Source": ["Paroxetine", "Paroxetine", "Paroxetine", "Paroxetine", "Paroxetine", "Paroxetine", "Paroxetine",
                "Paroxetine", "Paroxetine", "Paroxetine", "Paroxetine", "Paroxetine"]})


secstr_distorted_df = pd.concat([dssp_distorted_df_no_parox, paroxetine_distorted_df], ignore_index=True)

########
# PLOT
########

sns.set_style('darkgrid')

g = sns.catplot(data=secstr_distorted_df, x='Secondary structure', y='Percentage', hue='Variant', col='Source',
                kind='bar', legend_out=False, legend=False, height=4, aspect=0.9)

ax = g.facet_axis(0, 0)

for c in ax.containers:
    labels = [f'{(v.get_height()):.1f}' for v in c]
    ax.bar_label(c, labels=labels, label_type='edge', fontsize=14)

ax1 = g.facet_axis(0, 1)

for c in ax1.containers:
    labels = [f'{(v.get_height()):.1f}' for v in c]
    ax1.bar_label(c, labels=labels, label_type='edge', fontsize=14)

# ax2 = g.facet_axis(0, 2)

# for c in ax2.containers:
#     labels = [f'{(v.get_height()):.1f}' for v in c]
#     ax2.bar_label(c, labels=labels, label_type='edge', fontsize=14)

for label in (ax.get_xticklabels() + ax.get_yticklabels()):
    label.set_fontsize(16)

for label in ax.get_xticklabels():
    ax.set_xticklabels(["Helix", "Distorted\nHelix", "Strand", "Distorted\nStrand", "Turn", "Unordered"])

for label in (ax1.get_xticklabels() + ax1.get_yticklabels()):
    label.set_fontsize(16)

# for label in (ax2.get_xticklabels() + ax2.get_yticklabels()):
#     label.set_fontsize(16)

ax.set_xlabel('', size=14)
ax.set_ylabel('Percentage of amino acids', size=15)

# g.fig.text(0.69, 0.87, 'C', fontsize=22)
# g.fig.text(0.37, 0.87, 'B', fontsize=22)
# g.fig.text(0.048, 0.87, 'A', fontsize=22)


ax.set_title("No paroxetine", size=14)

ax1.set_xlabel('', size=2)
# ax1.set_xlabel('Secondary structure elements', size=16)
ax1.set_ylabel(None, size=14)

ax1.set_title("With paroxetine", size=14)

# ax2.set_xlabel('', size=14)
# ax2.set_ylabel(None, size=14)

# ax2.set_title("", size=14)

# leg = ax.legend(prop={"size": 16})
# leg1 = ax1.legend(prop={"size": 16}, loc='upper center')
# leg2 = ax2.legend(prop={"size": 16})

# mng = plt.get_current_fig_manager()
# mng.full_screen_toggle()

plt.show()

# g.figure.savefig("secstr_comparison_wt_k_paroxetine.png", dpi=320)












