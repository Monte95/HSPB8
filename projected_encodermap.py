import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
from remd_distributions import read_gyration_data

projected = np.load("projected.npy")

t, gyr_all_wt_500ns = read_gyration_data("gyrate_all_wt.xvg")

t, gyr_all_k_500ns = read_gyration_data("gyrate_all_k.xvg")

df_wt = pd.DataFrame({'Rg': gyr_all_wt_500ns})
df_k = pd.DataFrame({'Rg': gyr_all_k_500ns})

df_wt_k = pd.concat([df_wt, df_k], ignore_index=True)

df_wt_k['Variant'] = np.repeat(['wt', 'K141E'], int(len(gyr_all_wt_500ns)))

df_wt_k['Homology Model'] = np.repeat(['I-TASSER wt', 'MODELLER wt', 'ROSETTA wt', 'I-TASSER K141E',
                                       'MODELLER K141E', 'ROSETTA K141E'], [5001, 5000, 5000, 5001, 5000, 5000])

df_wt_k['x'] = projected[:, 0]
df_wt_k['y'] = projected[:, 1]

df_wt_k.loc[df_wt_k['Rg'] > 2.1, 'Conformation'] = 'Open'
df_wt_k.loc[df_wt_k['Rg'] <= 2.1, 'Conformation'] = 'Closed'

print([x+1 for x in df_k.loc[df_k['Rg'] <= 2.1].index.to_list()])

sns.set_style('darkgrid')
fig, ax = plt.subplots(1, figsize=(50, 20), gridspec_kw=dict(width_ratios=[4]))

# SCATTERPLOT

# sns.scatterplot(data=df_wt_k, x='x', y='y', hue='Conformation', s=10)
sns.scatterplot(data=df_wt_k, x='x', y='y', hue='Homology Model', s=10)

# ax.set_title("Open/Closed conformations", fontsize=20)
ax.set_title("Homology Model", fontsize=20)

ax.set_xlabel('X', fontsize=20)
ax.set_ylabel('Y', fontsize=20)

for label in (ax.get_xticklabels() + ax.get_yticklabels()):
    label.set_fontsize(14)

leg = ax.legend(prop={"size": 20})

# plt.show()

# fig.savefig('scatter_open_close_wt_k_emap.png', dpi=320)
fig.savefig('scatter_hm_wt_k_emap.png', dpi=320)
