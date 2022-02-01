from scipy import stats
import numpy as np
import matplotlib.pyplot as plt
from paroxetine_contacts import wt_k_paroxetine_contacts_frequency_df, wt_paroxetine_contacts_time_df, k_paroxetine_contacts_time_df

x_wt = np.array(wt_paroxetine_contacts_time_df['Residue'])
kde_wt = stats.gaussian_kde(x_wt)

x_k = np.array(k_paroxetine_contacts_time_df['Residue'])
kde_k = stats.gaussian_kde(x_k)

fig, ax = plt.subplots(figsize=(50, 40), nrows=2)
bins = np.arange(196)

N, bins, patches = ax[0].hist(x_wt, density=True, bins=bins, alpha=0.6)

for i in range(0, 96):
    patches[i].set_facecolor('blue')
for i in range(96, 170):
    patches[i].set_facecolor('yellow')
for i in range(170, 195):
    patches[i].set_facecolor('red')

ax[0].set(ylim=(0, 0.030))
ax[0].set_xlabel('', size=16)
ax[0].set_ylabel('Relative Frequency', size=16)
ax[0].set_title('wt', size=16)

xx = np.linspace(0, 196, 1000)
ax[0].plot(xx, kde_wt(xx), color='blue')

N, bins, patches = ax[1].hist(x_k, density=True, bins=bins, alpha=0.6)

for i in range(0, 96):
    patches[i].set_facecolor('blue')
for i in range(96, 170):
    patches[i].set_facecolor('yellow')
for i in range(170, 195):
    patches[i].set_facecolor('red')

ax[1].set_xlabel('Residue Number', size=16)
ax[1].set_ylabel('Relative Frequency', size=16)
ax[1].set_title('K141E', size=16)
ax[1].set(ylim=(0, 0.030))

ax[1].plot(xx, kde_k(xx), color='orange')

plt.show()

fig.savefig("paroxetine_contacts_residue_kde_wt_k.png", dpi=320)
