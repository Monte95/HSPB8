import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import pandas as pd
from scipy.signal import find_peaks


def read_gyration_data(document_name):
    x, y = [], []
    with open(document_name) as f:
        for line in f:
            cols = line.split()

            if len(cols) == 5:
                x.append(float(cols[0]))
                y.append(float(cols[1]))
    return np.array(x), np.array(y)


def read_dihedral_data(document_name):
    x, y = [], []
    with open(document_name) as f:
        for line in f:
            cols = line.split()

            if len(cols) == 2:
                x.append(float(cols[0]))
                y.append(abs(float(cols[1])))
    return x, y


def running_avg(y, num_ns, window_size=100):
    rmsd_series = pd.Series(y)
    windows = rmsd_series.rolling(window_size)
    rmsd_moving_averages = windows.mean()

    rmsd_moving_averages_list = rmsd_moving_averages.tolist()
    rmsd_moving_averages_without_nans = rmsd_moving_averages_list[window_size - 1:]

    x_running_avg = np.linspace(0, num_ns, len(rmsd_moving_averages_without_nans))

    return x_running_avg, rmsd_moving_averages_without_nans


def plot_gyration(num, x, y, labels, colors, ylim, figure_name):
    sns.set_style('darkgrid')

    fig = plt.figure()
    ax = fig.add_subplot(111)

    ax.set_title("Gyration radius", fontsize=20)
    ax.set_xlabel('Time (ns)', fontsize=20)
    ax.set_ylabel('Gyration radius (nm)', fontsize=20)

    plt.xticks(fontsize=20)
    plt.yticks(fontsize=20)
    plt.rc('legend', fontsize=15)

    plt.ylim(ylim[0], ylim[1])

    for i in range(0, num):
        ax.plot(x[i], y[i], color=colors[i], alpha=0.6, label=labels[i], linewidth=3)

    plt.hlines(2.1, 0, 500, linestyles='dashed')

    leg = ax.legend()

    plt.show()

    fig.savefig(figure_name)


def plot_kde_cumulative_seaborn(data1, data2, figname):
    df_wt = pd.DataFrame({'Rg': data1})
    df_k = pd.DataFrame({'Rg': data2})

    df_wt_k = pd.concat([df_wt, df_k], ignore_index=True)
    df_wt_k['Variant'] = np.repeat(['wt', 'K141E'], int(len(data1)))

    sns.set_style('darkgrid')

    fig, ax = plt.subplots(1, figsize=(50, 20), gridspec_kw=dict(width_ratios=[4]))

    sns.kdeplot(data=df_wt_k, x="Rg", hue="Variant", bw_adjust=0.8, linewidth=3, shade=True, common_norm=True,
                common_grid=True, ax=ax)

    plt.axvline([2.1], color='red', linestyle='--', linewidth=2)

    sns.set(font_scale=3)

    ax.set_title("Gyration radius", fontsize=30)
    ax.set_xlabel('Gyration radius (nm)', fontsize=30)
    ax.set_ylabel('Relative frequency', fontsize=30)

    plt.show()

    fig.savefig('{}.png'.format(figname), dpi=320)


def plot_kde_Rg_seaborn(data1, data2, data3, data4, figname):
    df_wt = pd.DataFrame({'Rg': data1})
    df_k = pd.DataFrame({'Rg': data2})
    df_wt_paroxetine = pd.DataFrame({'Rg': data3})
    df_k_paroxetine = pd.DataFrame({'Rg': data4})

    df_wt_k = pd.concat([df_wt, df_k], ignore_index=True)
    df_wt_k['Variant'] = np.repeat(['wt', 'K141E'], int(len(data1)))

    df_wt_k_paroxetine = pd.concat([df_wt_paroxetine, df_k_paroxetine], ignore_index=True)
    df_wt_k_paroxetine['Variant'] = np.repeat(['wt', 'K141E'], int(len(data3)))

    sns.set_style('darkgrid')

    fig, ax = plt.subplots(figsize=(50, 20), nrows=2)

    ax[0].axvline([2.1], color='red', linestyle='--', linewidth=2)

    sns.kdeplot(data=df_wt_k, x="Rg", hue="Variant", bw_adjust=0.8, linewidth=3, shade=True, common_norm=True,
                common_grid=True, ax=ax[0])

    ax[1].axvline([2.1], color='red', linestyle='--', linewidth=2)

    sns.kdeplot(data=df_wt_k_paroxetine, x="Rg", hue="Variant", bw_adjust=0.8, linewidth=3, shade=True, common_norm=True,
                common_grid=True, ax=ax[1])

    ax[0].set_title("No Paroxetine", fontsize=16)
    ax[0].set_xlabel('', fontsize=16)
    ax[0].set_ylabel('Relative frequency', fontsize=16)
    ax[0].set(xlim=(1.7, 3.5))

    ax[1].set_title("With Paroxetine", fontsize=16)
    ax[1].set_xlabel('Gyration Radius (nm)', fontsize=16)
    ax[1].set_ylabel('Relative frequency', fontsize=16)
    ax[1].set(xlim=(1.7, 3.5))

    plt.show()

    fig.savefig('{}.png'.format(figname), dpi=320)


def plot_kde_dihedrals_seaborn(data1, data2, gyr_wt, gyr_k, figname):
    df_wt = pd.DataFrame({'Dihedral 91-94': data1, 'Rg': gyr_wt})
    df_k = pd.DataFrame({'Dihedral 91-94': data2, 'Rg': gyr_k})

    df_wt_k = pd.concat([df_wt, df_k], ignore_index=True)
    df_wt_k['Variant'] = np.repeat(['wt', 'K141E'], int(len(data1)))

    df_wt_k.loc[df_wt_k['Rg'] > 2.1, 'Conformation'] = 'Open'
    df_wt_k.loc[df_wt_k['Rg'] <= 2.1, 'Conformation'] = 'Closed'

    sns.set_style('darkgrid')

    fig, ax = plt.subplots(1, figsize=(50, 20), gridspec_kw=dict(width_ratios=[4]))

    g = sns.kdeplot(data=df_wt_k, x='Dihedral 91-94', hue="Conformation", palette=['orange', 'blue'], bw_adjust=0.8,
                    linewidth=3, shade=True, clip=(0.0, 180.0), ax=ax)

    ax.set_title("Dihedral 91-94", fontsize=16)
    ax.set_xlabel('Dihedral angle', fontsize=16)
    ax.set_ylabel('Relative frequency', fontsize=16)

    ax.set_xticklabels(g.get_xticks(), size=16)
    ax.set_yticklabels(g.get_yticks(), size=16)

    plt.show()

    fig.savefig('{}.png'.format(figname), dpi=320)


t_wt_7_parox, gyr_wt_7_parox = read_gyration_data("gyrate_wt_7_paroxetine.xvg")
t_wt_8_parox, gyr_wt_8_parox = read_gyration_data("gyrate_wt_8_paroxetine.xvg")
t_wt_9_parox, gyr_wt_9_parox = read_gyration_data("gyrate_wt_9_paroxetine.xvg")
t_wt_10_parox, gyr_wt_10_parox = read_gyration_data("gyrate_wt_10_paroxetine.xvg")
# t_wt_11_parox, gyr_wt_11_parox = read_gyration_data("gyrate_wt_11_paroxetine.xvg")
t_wt_12_parox, gyr_wt_12_parox = read_gyration_data("gyrate_wt_12_paroxetine.xvg")
t_wt_13_parox, gyr_wt_13_parox = read_gyration_data("gyrate_wt_13_paroxetine.xvg")
t_wt_14_parox, gyr_wt_14_parox = read_gyration_data("gyrate_wt_14_paroxetine.xvg")
t_wt_15_parox, gyr_wt_15_parox = read_gyration_data("gyrate_wt_15_paroxetine.xvg")
t_wt_16_parox, gyr_wt_16_parox = read_gyration_data("gyrate_wt_16_paroxetine.xvg")
# t_wt_17_parox, gyr_wt_17_parox = read_gyration_data("gyrate_wt_17_paroxetine.xvg")
# t_wt_18_parox, gyr_wt_18_parox = read_gyration_data("gyrate_wt_18_paroxetine.xvg")
# t_wt_19_parox, gyr_wt_19_parox = read_gyration_data("gyrate_wt_19_paroxetine.xvg")

t_k_1_parox, gyr_k_1_parox = read_gyration_data("gyrate_k_1_paroxetine.xvg")
t_k_2_parox, gyr_k_2_parox = read_gyration_data("gyrate_k_2_paroxetine.xvg")
t_k_3_parox, gyr_k_3_parox = read_gyration_data("gyrate_k_3_paroxetine.xvg")
t_k_4_parox, gyr_k_4_parox = read_gyration_data("gyrate_k_4_paroxetine.xvg")
t_k_5_parox, gyr_k_5_parox = read_gyration_data("gyrate_k_5_paroxetine.xvg")
t_k_6_parox, gyr_k_6_parox = read_gyration_data("gyrate_k_6_paroxetine.xvg")
t_k_20_parox, gyr_k_20_parox = read_gyration_data("gyrate_k_20_paroxetine.xvg")
t_k_21_parox, gyr_k_21_parox = read_gyration_data("gyrate_k_21_paroxetine.xvg")
t_k_22_parox, gyr_k_22_parox = read_gyration_data("gyrate_k_22_paroxetine.xvg")

t_wt_parox, gyr_wt_parox = read_gyration_data("gyrate_wt_paroxetine.xvg")
t_k_parox, gyr_k_parox = read_gyration_data("gyrate_k_paroxetine.xvg")

t_wt, gyr_all_wt = read_gyration_data("gyrate_all_wt.xvg")
t_k, gyr_all_k = read_gyration_data("gyrate_all_k.xvg")

# GYRATION ALL WT + K
t_it_wt, gyr_it_wt = read_gyration_data("../gyrate_it_wt_remd_500ns.xvg")
t_m_wt, gyr_m_wt = read_gyration_data("../gyrate_m_wt_remd_500ns.xvg")
t_trros_wt, gyr_trros_wt = read_gyration_data("../gyrate_trros_wt_remd_500ns.xvg")

t_it_k, gyr_it_k = read_gyration_data("../gyrate_it_k_remd_500ns.xvg")
t_m_k, gyr_m_k = read_gyration_data("../gyrate_m_k_remd_500ns.xvg")
t_trros_k, gyr_trros_k = read_gyration_data("../gyrate_trros_k_remd_500ns.xvg")

# FRAMES OPEN/CLOSED
df_wt = pd.DataFrame({'Rg': gyr_wt_parox})
df_k = pd.DataFrame({'Rg': gyr_k_parox})

# EXTRACT OPEN/CLOSED FRAMES
#print([x+1 for x in df_wt.loc[df_wt['Rg'] <= 2.1].index.to_list()])
print(len(df_wt.loc[df_wt['Rg'] > 2.7])/(len(df_wt['Rg']))*100)
print(len(df_k.loc[df_k['Rg'] <= 2.1])/(len(df_k['Rg']))*100)


# RUNNING AVERAGES
t_wt_7_parox_runavg, gyr_wt_7_parox_runavg = running_avg(gyr_wt_7_parox, len(gyr_wt_7_parox) / 10)
t_wt_8_parox_runavg, gyr_wt_8_parox_runavg = running_avg(gyr_wt_8_parox, len(gyr_wt_8_parox) / 10)
t_wt_9_parox_runavg, gyr_wt_9_parox_runavg = running_avg(gyr_wt_9_parox, len(gyr_wt_9_parox) / 10)
t_wt_10_parox_runavg, gyr_wt_10_parox_runavg = running_avg(gyr_wt_10_parox, len(gyr_wt_10_parox) / 10)
# t_wt_11_parox_runavg, gyr_wt_11_parox_runavg = running_avg(gyr_wt_11_parox, len(gyr_wt_11_parox) / 10)
t_wt_12_parox_runavg, gyr_wt_12_parox_runavg = running_avg(gyr_wt_12_parox, len(gyr_wt_12_parox) / 10)
t_wt_13_parox_runavg, gyr_wt_13_parox_runavg = running_avg(gyr_wt_13_parox, len(gyr_wt_13_parox) / 10)
t_wt_14_parox_runavg, gyr_wt_14_parox_runavg = running_avg(gyr_wt_14_parox, len(gyr_wt_14_parox) / 10)
t_wt_15_parox_runavg, gyr_wt_15_parox_runavg = running_avg(gyr_wt_15_parox, len(gyr_wt_15_parox) / 10)
t_wt_16_parox_runavg, gyr_wt_16_parox_runavg = running_avg(gyr_wt_16_parox, len(gyr_wt_16_parox) / 10)
# t_wt_17_parox_runavg, gyr_wt_17_parox_runavg = running_avg(gyr_wt_17_parox, len(gyr_wt_17_parox) / 10)
# t_wt_18_parox_runavg, gyr_wt_18_parox_runavg = running_avg(gyr_wt_18_parox, len(gyr_wt_18_parox) / 10)
# t_wt_19_parox_runavg, gyr_wt_19_parox_runavg = running_avg(gyr_wt_19_parox, len(gyr_wt_19_parox) / 10)

t_k_1_parox_runavg, gyr_k_1_parox_runavg = running_avg(gyr_k_1_parox, len(gyr_k_4_parox) / 10)
t_k_2_parox_runavg, gyr_k_2_parox_runavg = running_avg(gyr_k_2_parox, len(gyr_k_2_parox) / 10)
t_k_3_parox_runavg, gyr_k_3_parox_runavg = running_avg(gyr_k_3_parox, len(gyr_k_3_parox) / 10)
t_k_4_parox_runavg, gyr_k_4_parox_runavg = running_avg(gyr_k_4_parox, len(gyr_k_4_parox) / 10)
t_k_5_parox_runavg, gyr_k_5_parox_runavg = running_avg(gyr_k_5_parox, len(gyr_k_5_parox) / 10)
t_k_6_parox_runavg, gyr_k_6_parox_runavg = running_avg(gyr_k_6_parox, len(gyr_k_6_parox) / 10)
t_k_20_parox_runavg, gyr_k_20_parox_runavg = running_avg(gyr_k_20_parox, len(gyr_k_20_parox) / 10)
t_k_21_parox_runavg, gyr_k_21_parox_runavg = running_avg(gyr_k_21_parox, len(gyr_k_21_parox) / 10)
t_k_22_parox_runavg, gyr_k_22_parox_runavg = running_avg(gyr_k_22_parox, len(gyr_k_22_parox) / 10)

gyration_wt = [gyr_it_wt, gyr_m_wt, gyr_trros_wt]

gyration_wt_paroxetine = [gyr_wt_7_parox, gyr_wt_8_parox, gyr_wt_9_parox, gyr_wt_10_parox,
                          gyr_wt_12_parox, gyr_wt_13_parox, gyr_wt_14_parox,
                          gyr_wt_15_parox, gyr_wt_16_parox]

time_runavg_wt_paroxetine = [t_wt_7_parox_runavg, t_wt_8_parox_runavg, t_wt_9_parox_runavg, t_wt_10_parox_runavg,
                             t_wt_12_parox_runavg, t_wt_13_parox_runavg, t_wt_14_parox_runavg,
                             t_wt_15_parox_runavg, t_wt_16_parox_runavg]

gyration_runavg_wt_paroxetine = [gyr_wt_7_parox_runavg, gyr_wt_8_parox_runavg, gyr_wt_9_parox_runavg,
                                 gyr_wt_10_parox_runavg,
                                 gyr_wt_12_parox_runavg, gyr_wt_13_parox_runavg, gyr_wt_14_parox_runavg,
                                 gyr_wt_15_parox_runavg, gyr_wt_16_parox_runavg]

gyration_k = [gyr_it_k, gyr_m_k, gyr_trros_k]

gyration_k_paroxetine = [gyr_k_1_parox, gyr_k_2_parox, gyr_k_3_parox, gyr_k_4_parox,
                         gyr_k_5_parox, gyr_k_6_parox, gyr_k_20_parox, gyr_k_21_parox,
                         gyr_k_22_parox]

time_runavg_k_paroxetine = [t_k_1_parox_runavg, t_k_2_parox_runavg, t_k_3_parox_runavg, t_k_4_parox_runavg,
                            t_k_5_parox_runavg, t_k_6_parox_runavg, t_k_20_parox_runavg, t_k_21_parox_runavg,
                            t_k_22_parox_runavg]

gyration_runavg_k_paroxetine = [gyr_k_1_parox_runavg, gyr_k_2_parox_runavg, gyr_k_3_parox_runavg, gyr_k_4_parox_runavg,
                                gyr_k_5_parox_runavg, gyr_k_6_parox_runavg, gyr_k_20_parox_runavg,
                                gyr_k_21_parox_runavg,
                                gyr_k_22_parox_runavg]

labels_gyration_comparison_wt = ['wt 7', 'wt 8', 'wt 9', 'wt 10', 'wt 12', 'wt 13', 'wt 14', 'wt 15', 'wt 16']
labels_gyration_comparison_k = ['K141E 1', 'K141E 2', 'K141E 3', 'K141E 4', 'K141E 5', 'K141E 6', 'K141E 20',
                                'K141E 21', 'K141E 22']
colors = ['blue', 'orange', 'green', 'red', 'cyan', 'purple', 'magenta', 'black', 'yellow']

###########
# DIHEDRALS
time_wt_paroxetine, dihedral_91_94_wt_paroxetine = read_dihedral_data("dihedrals_91_94_wt_paroxetine.dat")
time_k_paroxetine, dihedral_91_94_k_paroxetine = read_dihedral_data("dihedrals_91_94_k_paroxetine.dat")

plot_kde_dihedrals_seaborn(dihedral_91_94_wt_paroxetine[0:-1], dihedral_91_94_k_paroxetine[0:-1], gyr_wt_parox, gyr_k_parox,
                               figname='dihedrals_91_94_kde_wt_k_paroxetine')

# df_wt = pd.DataFrame({'Dihedral 91-94': dihedral_91_94_wt_paroxetine[0:-1], 'Rg': gyr_wt_parox})
# df_k = pd.DataFrame({'Dihedral 91-94': dihedral_91_94_k_paroxetine[0:-1], 'Rg': gyr_k_parox})
#
#
# df_wt_k = pd.concat([df_wt, df_k], ignore_index=True)
# df_wt_k['Variant'] = np.repeat(['wt', 'K141E'], int(len(dihedral_91_94_wt_paroxetine[:-1])))
#
# df_wt_k.loc[df_wt_k['Rg'] <= 2.1, 'Conformation'] = 'Closed'
# df_wt_k.loc[df_wt_k['Rg'] > 2.1, 'Conformation'] = 'Open'
#
# sns.set_style('darkgrid')
#
# fig, ax = plt.subplots(1, figsize=(50, 20), gridspec_kw=dict(width_ratios=[4]))
#
# g = sns.kdeplot(data=df_wt_k, x='Dihedral 91-94', hue="Conformation", bw_adjust=0.8, clip=(0.0, 180.0), ax=ax)
#
# plt.show()
#
# x = ax.lines[0].get_xdata()  # Get the x data of the distribution
# y = ax.lines[0].get_ydata()  # Get the y data of the distribution
# peaks, _ = np.array(find_peaks(y))
#
# for p in range(len(peaks)):
#     if y[peaks[p]] >= 0.001:
#         print(x[peaks[p]])
#         print(y[peaks[p]])

##############################################
# COUNT OPEN/CLOSED AND CLOSED/OPEN TRANSITIONS
count_oc_wt = 0
count_co_wt = 0

count_oc_k = 0
count_co_k = 0

ns_wt = len(gyration_wt[0])+len(gyration_wt[1])+len(gyration_wt[2])

num_100ns_wt = np.sum([len(gyr_wt) for gyr_wt in gyration_wt]) * 500 / (5001 * 100)
num_100ns_k = np.sum([len(gyr_k) for gyr_k in gyration_k]) * 500 / (5001 * 100)

num_100ns_wt_paroxetine = np.sum([len(gyr_wt) for gyr_wt in gyration_wt_paroxetine]) * 500 / (5001 * 100)
num_100ns_k_paroxetine = np.sum([len(gyr_k) for gyr_k in gyration_k_paroxetine]) * 500 / (5001 * 100)

for gyr_wt in gyration_wt:
    for j in range(1, len(gyr_wt)):
        if gyr_wt[j] < 2.1 and gyr_wt[j - 1] > 2.1:
            count_oc_wt += 1

        if gyr_wt[j] > 2.1 and gyr_wt[j - 1] < 2.1:
            count_co_wt += 1

for gyr_k in gyration_k:
    for j in range(1, len(gyr_k)):
        if gyr_k[j] < 2.1 and gyr_k[j - 1] > 2.1:
            count_oc_k += 1

        if gyr_k[j] > 2.1 and gyr_k[j - 1] < 2.1:
            count_co_k += 1

count_oc_wt_100ns = count_oc_wt / num_100ns_wt
count_co_wt_100ns = count_co_wt / num_100ns_wt

count_oc_k_100ns = count_oc_k / num_100ns_k
count_co_k_100ns = count_co_k / num_100ns_k

count_oc_wt_paroxetine = 0
count_co_wt_paroxetine = 0

count_oc_k_paroxetine = 0
count_co_k_paroxetine = 0

for gyr_wt in gyration_wt_paroxetine:
    for j in range(1, len(gyr_wt)):
        if gyr_wt[j] < 2.1 and gyr_wt[j - 1] > 2.1:
            count_oc_wt_paroxetine += 1

        if gyr_wt[j] > 2.1 and gyr_wt[j - 1] < 2.1:
            count_co_wt_paroxetine += 1

for gyr_k in gyration_k_paroxetine:
    for j in range(1, len(gyr_k)):
        if gyr_k[j] < 2.1 and gyr_k[j - 1] > 2.1:
            count_oc_k_paroxetine += 1

        if gyr_k[j] > 2.1 and gyr_k[j - 1] < 2.1:
            count_co_k_paroxetine += 1

count_oc_wt_paroxetine_100ns = count_oc_wt_paroxetine / num_100ns_wt_paroxetine
count_co_wt_paroxetine_100ns = count_co_wt_paroxetine / num_100ns_wt_paroxetine

count_oc_k_paroxetine_100ns = count_oc_k_paroxetine / num_100ns_k_paroxetine
count_co_k_paroxetine_100ns = count_co_k_paroxetine / num_100ns_k_paroxetine

open_closed_transitions_df = pd.DataFrame({'Transition': ['Open -> Closed', 'Closed -> Open',
                                                          'Open -> Closed', 'Closed -> Open',
                                                          'Open -> Closed', 'Closed -> Open',
                                                          'Open -> Closed', 'Closed -> Open'],
                                           'Paroxetine': ['No paroxetine', 'No paroxetine',
                                                          'No paroxetine', 'No paroxetine',
                                                          'Paroxetine', 'Paroxetine',
                                                          'Paroxetine', 'Paroxetine'],
                                           'Count': [count_oc_wt_100ns, count_co_wt_100ns,
                                                     count_oc_k_100ns, count_co_k_100ns,
                                                     count_oc_wt_paroxetine_100ns, count_co_wt_paroxetine_100ns,
                                                     count_oc_k_paroxetine_100ns, count_co_k_paroxetine_100ns],
                                           'Variant': ['wt', 'wt', 'K141E', 'K141E', 'wt', 'wt', 'K141E', 'K141E']})

# PLOT OPEN/CLOSED TRANSITIONS
sns.set_style('darkgrid')

g = sns.catplot(data=open_closed_transitions_df, x='Transition', y='Count', hue='Variant', col='Paroxetine',
                kind='bar', legend_out=False, legend=False)

g.set(ylim=(0, 25))

ax = g.facet_axis(0, 0)

for c in ax.containers:
    labels = [f'{(v.get_height()):.0f}' for v in c]
    ax.bar_label(c, labels=labels, label_type='edge', fontsize=16)

for label in (ax.get_xticklabels() + ax.get_yticklabels()):
    label.set_fontsize(16)

ax.set_xlabel('', size=16)
ax.set_ylabel('Number of transitions in 100 ns', size=16)

ax.set_title("No Paroxetine", size=16)

leg = ax.legend(prop={"size": 18}, loc='upper right')

ax2 = g.facet_axis(0, 1)

for c in ax2.containers:
    labels = [f'{(v.get_height()):.0f}' for v in c]
    ax2.bar_label(c, labels=labels, label_type='edge', fontsize=16)

for label in (ax2.get_xticklabels() + ax2.get_yticklabels()):
    label.set_fontsize(16)

ax2.set_xlabel('', size=16)
ax2.set_ylabel('Number of transitions in 100 ns', size=16)

ax2.set_title("With Paroxetine", size=16)

leg = ax2.legend(prop={"size": 18}, loc='upper right')

plt.close()
# plt.show()

# g.figure.savefig("open_closed_transitions_wt_k_paroxetine.png", dpi=320)

# plot_gyration(9, time_runavg_k, gyration_runavg_k, labels_gyration_comparison_k,
#               colors, [1.7, 3.25], 'gyration_k_parox_runavg')


# KDE
# plot_kde_cumulative_seaborn(gyr_wt_parox, gyr_k_parox,
#                             figname='gyration_wt_k_parox')

# plt.show()

############################
# PLOT RG TREMD + PAROXETINE
# plot_kde_Rg_seaborn(gyr_all_wt, gyr_all_k, gyr_wt_parox, gyr_k_parox, figname='gyration_wt_k_parox')

