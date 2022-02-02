import mdtraj as md
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import time

NUM_RESIDUES = 196

ref_structure_wt = "frame1_wt_paroxetine_open.gro"
trajectory_wt = "wt_paroxetine.xtc"

ref_structure_k = "frame1_k_paroxetine_open.gro"
trajectory_k = "k_paroxetine.xtc"

t_wt = md.load(trajectory_wt, top=ref_structure_wt)
t_k = md.load(trajectory_k, top=ref_structure_k)

start_time = time.time()

dssp_wt = np.transpose(md.compute_dssp(t_wt, simplified=False))
dssp_k = np.transpose(md.compute_dssp(t_k, simplified=False))

print("--- %s seconds ---" % (time.time() - start_time))

dssp_wt_df = pd.DataFrame(dssp_wt)
dssp_k_df = pd.DataFrame(dssp_k)

np.set_printoptions(threshold=np.inf)

secstr_wt_tot_dssp = (dssp_wt_df.apply(pd.value_counts).sum(axis=1) / (NUM_RESIDUES * len(t_wt))) * 100
secstr_k_tot_dssp = (dssp_k_df.apply(pd.value_counts).sum(axis=1) / (NUM_RESIDUES * len(t_k))) * 100

print(secstr_wt_tot_dssp)
print(secstr_k_tot_dssp)

secstr_wt_residues_dssp = (dssp_wt_df.apply(pd.Series.value_counts, axis=1).fillna(0) / len(t_wt)) * 100
secstr_k_residues_dssp = (dssp_k_df.apply(pd.Series.value_counts, axis=1).fillna(0) / len(t_k)) * 100

#######
# PLOT
#######

sns.set_style('darkgrid')
fig, ax = plt.subplots(nrows=2, ncols=1, figsize=(60, 30), sharex=True)

# ax[1].xaxis.set_tick_params(labelsize=10)

columns = secstr_wt_residues_dssp.columns
colors_3 = ['#3a86ff', '#ffbe0b', '#8338ec']
colors_8 = ['#0d3b66', '#f48c06', '#ffbe0b', '#7209b7', '#8338ec', '#9d4edd', '#38a3a5', '#00bbf9']

a = secstr_wt_residues_dssp[columns].plot(kind="bar", stacked=True, color=colors_8, ax=ax[0])

ax[0].set_title("wt", fontsize=20)
ax[0].set_ylabel('Percentage', fontsize=20)

ax[0].xaxis.set_major_locator(ticker.AutoLocator())
ax[0].xaxis.set_major_formatter(ticker.ScalarFormatter(useMathText=True))

b = secstr_k_residues_dssp[columns].plot(kind="bar", stacked=True, color=colors_8, ax=ax[1])

plt.xticks(np.arange(1, 197, 5))

ax[1].set_title("K141E", fontsize=20)
ax[1].set_xlabel('Residue', fontsize=20)
ax[1].set_ylabel('Percentage', fontsize=20)

leg = ax[1].legend(loc='lower right')
# plt.show()

#fig.savefig('dssp_per_residue_8cat_wt_k.png', dpi=320)

# CD DATA
#           HELIX                           STRAND                          TURNS     UNORDERED
#           REGULAR   DISTORTED   TOTAL     REGULAR     DISTORTED   TOTAL
# wt        0.012     0.040       0.052     0.251       0.114       0.365   0.208     0.374
# K141E     0.012     0.045       0.057     0.228       0.092       0.320   0.165     0.458

# DSSP
# 3-CATEGORY ASSIGNMENT
#           HELIX    STRAND   COIL
# wt        0.098    0.202    0.699
# K141E     0.082    0.197    0.721

# DSSP
# 8-CATEGORY ASSIGNMENT
#          HELIX                            STRAND                      UNORDERED
#          ALPHA   3_10    PI      TOTAL    EXTENDED  BRIDGE   TOTAL    TURN    BEND    IRREGULAR   TOTAL
# wt       0.067   0.031   0.000   0.098    0.184     0.019    0.203    0.131   0.197   0.371       0.699
# K141E    0.048   0.033   0.000   0.081    0.183     0.014    0.197    0.141   0.187   0.393       0.721

# DSSP PAROXETINE
# 8-CATEGORY ASSIGNMENT
#          HELIX                            STRAND                      UNORDERED
#          ALPHA   3_10    PI      TOTAL    EXTENDED  BRIDGE   TOTAL    TURN    BEND    IRREGULAR   TOTAL
# wt       0.060   0.039   0.000   0.099    0.181     0.019    0.200    0.137   0.190   0.373       0.700
# K141E    0.054   0.036   0.000   0.090    0.189     0.016    0.205    0.137   0.175   0.391       0.703

