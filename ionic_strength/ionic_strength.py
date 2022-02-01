import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt


def read_gyration_data(document_name):
    x, y = [], []
    with open(document_name) as f:
        for line in f:
            cols = line.split()

            if len(cols) == 5:
                x.append(float(cols[0]))
                y.append(float(cols[1]))
    return np.array(x), np.array(y)


def plot_kde_Rg_seaborn(data1, data2, data3, data4, figname):
    df_wt = pd.DataFrame({'Rg': data1})
    df_k = pd.DataFrame({'Rg': data2})
    df_wt_50mM = pd.DataFrame({'Rg': data3})
    df_k_50mM = pd.DataFrame({'Rg': data4})

    df_wt_k = pd.concat([df_wt, df_k], ignore_index=True)
    df_wt_k['Variant'] = np.repeat(['wt', 'K141E'], int(len(data1)))

    df_wt_k_50mM = pd.concat([df_wt_50mM, df_k_50mM], ignore_index=True)
    df_wt_k_50mM['Variant'] = np.repeat(['wt', 'K141E'], int(len(data3)))

    sns.set_style('darkgrid')

    fig, ax = plt.subplots(figsize=(50, 20), nrows=2)

    ax[0].axvline([2.1], color='red', linestyle='--', linewidth=2)

    sns.kdeplot(data=df_wt_k, x="Rg", hue="Variant", bw_adjust=0.8, linewidth=3, shade=True, common_norm=True,
                common_grid=True, ax=ax[0])

    ax[1].axvline([2.1], color='red', linestyle='--', linewidth=2)

    sns.kdeplot(data=df_wt_k_50mM, x="Rg", hue="Variant", bw_adjust=0.8, linewidth=3, shade=True, common_norm=True,
                common_grid=True, ax=ax[1])

    ax[0].set_title("0 mM", fontsize=16)
    ax[0].set_xlabel('', fontsize=16)
    ax[0].set_ylabel('Relative frequency', fontsize=16)
    ax[0].set(xlim=(1.7, 3.5))

    ax[1].set_title("50 mM", fontsize=16)
    ax[1].set_xlabel('Gyration Radius (nm)', fontsize=16)
    ax[1].set_ylabel('Relative frequency', fontsize=16)
    ax[1].set(xlim=(1.7, 3.5))

    plt.show()

    fig.savefig('{}.png'.format(figname), dpi=320)

# IMPORT DATA
t_wt_8_50mM, gyr_wt_8_50mM = read_gyration_data("gyrate_wt_8_50mM.xvg")
# t_k_5_50mM, gyr_k_5_50mM = read_gyration_data("gyrate_k_5_50mM.xvg")

t_all_wt, gyr_all_wt = read_gyration_data("gyrate_all_wt.xvg")
t_all_k, gyr_all_k = read_gyration_data("gyrate_all_k.xvg")

# plot_kde_Rg_seaborn(gyr_all_wt, gyr_all_k, gyr_wt_8_50mM, gyr_k_5_50mM, figname='gyrate_wt_k_50mM')

df_wt = pd.DataFrame({'Rg': gyr_wt_8_50mM})

# PLOT KDE GYRATION RADIUS
sns.set_style('darkgrid')

fig, ax = plt.subplots(1, figsize=(50, 20), gridspec_kw=dict(width_ratios=[4]))

sns.kdeplot(data=df_wt, x="Rg", bw_adjust=0.8, linewidth=3, shade=True, common_norm=True,
                common_grid=True, ax=ax)

plt.axvline([2.1], color='red', linestyle='--', linewidth=2)
ax.set_title("", fontsize=30)
ax.set_xlabel('Gyration radius (nm)', fontsize=18)
ax.set_ylabel('Relative frequency', fontsize=18)

for label in (ax.get_xticklabels() + ax.get_yticklabels()):
    label.set_fontsize(16)

plt.show()
plt.close()
# fig.savefig('gyrate_wt_50mM.png', dpi=320)
