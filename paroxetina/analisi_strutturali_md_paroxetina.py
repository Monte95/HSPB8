import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns


def read_gyration_data(document_name):
    x, y = [], []
    with open(document_name) as f:
        for line in f:
            cols = line.split()

            if len(cols) == 5:
                x.append(float(cols[0]))
                y.append(float(cols[1]))
    return np.array(x), np.array(y)


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

def plot_kde_gyration(data1, data2, data3, data4, figname):
    df_wt = pd.DataFrame({'Rg': data1})
    df_k = pd.DataFrame({'Rg': data2})
    df_wt_paroxetine = pd.DataFrame({'Rg': data3})
    df_k_paroxetine = pd.DataFrame({'Rg': data4})

    df_wt_k = pd.concat([df_wt, df_k], ignore_index=True)
    df_wt_k['Variant'] = np.repeat(['wt', 'K141E'], int(len(data1)))

    df_wt_k_paroxetine = pd.concat([df_wt_paroxetine, df_k_paroxetine], ignore_index=True)
    df_wt_k_paroxetine['Variant'] = np.repeat(['wt', 'K141E'], [int(len(data3)), int(len(data4))])

    sns.set_style('darkgrid')

    fig, ax = plt.subplots(figsize=(50, 20), nrows=2)

    ax[0].axvline([2.1], color='red', linestyle='--', linewidth=2)

    a = sns.kdeplot(data=df_wt_k, x="Rg", hue="Variant", bw_adjust=0.8, linewidth=3, shade=True, common_norm=True,
                    common_grid=True, ax=ax[0], legend=False)

    ax[1].axvline([2.1], color='red', linestyle='--', linewidth=2)

    b = sns.kdeplot(data=df_wt_k_paroxetine, x="Rg", hue="Variant", bw_adjust=0.8, linewidth=3, shade=True,
                    common_norm=True,
                    common_grid=True, ax=ax[1], legend=False)

    ax[0].set_title("No Paroxetine", fontsize=16)
    ax[0].set_xlabel('', fontsize=16)
    ax[0].set_ylabel('Relative frequency', fontsize=20)
    ax[0].set(xlim=(1.7, 3.5))

    ax[1].set_title("With Paroxetine", fontsize=16)
    ax[1].set_xlabel('Gyration Radius (nm)', fontsize=20)
    ax[1].set_ylabel('Relative frequency', fontsize=20)
    ax[1].set(xlim=(1.7, 3.5))

    for label in (ax[0].get_xticklabels() + ax[0].get_yticklabels()):
        label.set_fontsize(16)

    for label in (ax[1].get_xticklabels() + ax[1].get_yticklabels()):
        label.set_fontsize(16)

    plt.subplots_adjust(top=0.96, bottom=0.08, left=0.049, right=0.992, hspace=0.175, wspace=0.2)

    plt.show()

    fig.savefig('../images/{}.png'.format(figname), dpi=320)


# READ GYRATION DATA
t_wt_7_parox, gyr_wt_7_parox = read_gyration_data("../paroxetine/gyrate_wt_7_paroxetine.xvg")
t_wt_8_parox, gyr_wt_8_parox = read_gyration_data("../paroxetine/gyrate_wt_8_paroxetine.xvg")
t_wt_9_parox, gyr_wt_9_parox = read_gyration_data("../paroxetine/gyrate_wt_9_paroxetine.xvg")
t_wt_10_parox, gyr_wt_10_parox = read_gyration_data("../paroxetine/gyrate_wt_10_paroxetine.xvg")
t_wt_11_parox, gyr_wt_11_parox = read_gyration_data("../paroxetine/gyrate_wt_11_paroxetine.xvg")
t_wt_12_parox, gyr_wt_12_parox = read_gyration_data("../paroxetine/gyrate_wt_12_paroxetine.xvg")
t_wt_13_parox, gyr_wt_13_parox = read_gyration_data("../paroxetine/gyrate_wt_13_paroxetine.xvg")
t_wt_14_parox, gyr_wt_14_parox = read_gyration_data("../paroxetine/gyrate_wt_14_paroxetine.xvg")
t_wt_15_parox, gyr_wt_15_parox = read_gyration_data("../paroxetine/gyrate_wt_15_paroxetine.xvg")
t_wt_16_parox, gyr_wt_16_parox = read_gyration_data("../paroxetine/gyrate_wt_16_paroxetine.xvg")
t_wt_17_parox, gyr_wt_17_parox = read_gyration_data("../paroxetine/gyrate_wt_17_paroxetine.xvg")
t_wt_18_parox, gyr_wt_18_parox = read_gyration_data("../paroxetine/gyrate_wt_18_paroxetine.xvg")
t_wt_19_parox, gyr_wt_19_parox = read_gyration_data("../paroxetine/gyrate_wt_19_paroxetine.xvg")

t_k_1_parox, gyr_k_1_parox = read_gyration_data("../paroxetine/gyrate_k_1_paroxetine.xvg")
t_k_2_parox, gyr_k_2_parox = read_gyration_data("../paroxetine/gyrate_k_2_paroxetine.xvg")
t_k_3_parox, gyr_k_3_parox = read_gyration_data("../paroxetine/gyrate_k_3_paroxetine.xvg")
t_k_4_parox, gyr_k_4_parox = read_gyration_data("../paroxetine/gyrate_k_4_paroxetine.xvg")
t_k_5_parox, gyr_k_5_parox = read_gyration_data("../paroxetine/gyrate_k_5_paroxetine.xvg")
t_k_6_parox, gyr_k_6_parox = read_gyration_data("../paroxetine/gyrate_k_6_paroxetine.xvg")
t_k_20_parox, gyr_k_20_parox = read_gyration_data("../paroxetine/gyrate_k_20_paroxetine.xvg")
t_k_21_parox, gyr_k_21_parox = read_gyration_data("../paroxetine/gyrate_k_21_paroxetine.xvg")
t_k_22_parox, gyr_k_22_parox = read_gyration_data("../paroxetine/gyrate_k_22_paroxetine.xvg")

t_wt_parox, gyr_wt_parox = read_gyration_data("../paroxetine/gyrate_wt_paroxetine_all.xvg")
t_k_parox, gyr_k_parox = read_gyration_data("../paroxetine/gyrate_k_paroxetine.xvg")

t_wt, gyr_all_wt = read_gyration_data("../tremd/gyrate_all_wt.xvg")
t_k, gyr_all_k = read_gyration_data("../tremd/gyrate_all_k.xvg")

# RUNNING AVERAGES
t_wt_7_parox_runavg, gyr_wt_7_parox_runavg = running_avg(gyr_wt_7_parox, len(gyr_wt_7_parox) / 10)
t_wt_8_parox_runavg, gyr_wt_8_parox_runavg = running_avg(gyr_wt_8_parox, len(gyr_wt_8_parox) / 10)
t_wt_9_parox_runavg, gyr_wt_9_parox_runavg = running_avg(gyr_wt_9_parox, len(gyr_wt_9_parox) / 10)
t_wt_10_parox_runavg, gyr_wt_10_parox_runavg = running_avg(gyr_wt_10_parox, len(gyr_wt_10_parox) / 10)
t_wt_11_parox_runavg, gyr_wt_11_parox_runavg = running_avg(gyr_wt_11_parox, len(gyr_wt_11_parox) / 10)
t_wt_12_parox_runavg, gyr_wt_12_parox_runavg = running_avg(gyr_wt_12_parox, len(gyr_wt_12_parox) / 10)
t_wt_13_parox_runavg, gyr_wt_13_parox_runavg = running_avg(gyr_wt_13_parox, len(gyr_wt_13_parox) / 10)
t_wt_14_parox_runavg, gyr_wt_14_parox_runavg = running_avg(gyr_wt_14_parox, len(gyr_wt_14_parox) / 10)
t_wt_15_parox_runavg, gyr_wt_15_parox_runavg = running_avg(gyr_wt_15_parox, len(gyr_wt_15_parox) / 10)
t_wt_16_parox_runavg, gyr_wt_16_parox_runavg = running_avg(gyr_wt_16_parox, len(gyr_wt_16_parox) / 10)
t_wt_17_parox_runavg, gyr_wt_17_parox_runavg = running_avg(gyr_wt_17_parox, len(gyr_wt_17_parox) / 10)
t_wt_18_parox_runavg, gyr_wt_18_parox_runavg = running_avg(gyr_wt_18_parox, len(gyr_wt_18_parox) / 10)
t_wt_19_parox_runavg, gyr_wt_19_parox_runavg = running_avg(gyr_wt_19_parox, len(gyr_wt_19_parox) / 10)

t_k_1_parox_runavg, gyr_k_1_parox_runavg = running_avg(gyr_k_1_parox, len(gyr_k_4_parox) / 10)
t_k_2_parox_runavg, gyr_k_2_parox_runavg = running_avg(gyr_k_2_parox, len(gyr_k_2_parox) / 10)
t_k_3_parox_runavg, gyr_k_3_parox_runavg = running_avg(gyr_k_3_parox, len(gyr_k_3_parox) / 10)
t_k_4_parox_runavg, gyr_k_4_parox_runavg = running_avg(gyr_k_4_parox, len(gyr_k_4_parox) / 10)
t_k_5_parox_runavg, gyr_k_5_parox_runavg = running_avg(gyr_k_5_parox, len(gyr_k_5_parox) / 10)
t_k_6_parox_runavg, gyr_k_6_parox_runavg = running_avg(gyr_k_6_parox, len(gyr_k_6_parox) / 10)
t_k_20_parox_runavg, gyr_k_20_parox_runavg = running_avg(gyr_k_20_parox, len(gyr_k_20_parox) / 10)
t_k_21_parox_runavg, gyr_k_21_parox_runavg = running_avg(gyr_k_21_parox, len(gyr_k_21_parox) / 10)
t_k_22_parox_runavg, gyr_k_22_parox_runavg = running_avg(gyr_k_22_parox, len(gyr_k_22_parox) / 10)

gyration_wt_paroxetine = [gyr_wt_7_parox, gyr_wt_8_parox, gyr_wt_9_parox, gyr_wt_10_parox, gyr_wt_11_parox,
                          gyr_wt_12_parox, gyr_wt_13_parox, gyr_wt_14_parox, gyr_wt_15_parox, gyr_wt_16_parox,
                          gyr_wt_17_parox, gyr_wt_18_parox, gyr_wt_19_parox]

time_runavg_wt_paroxetine = [t_wt_7_parox_runavg, t_wt_8_parox_runavg, t_wt_9_parox_runavg, t_wt_10_parox_runavg,
                             t_wt_11_parox_runavg, t_wt_12_parox_runavg, t_wt_13_parox_runavg, t_wt_14_parox_runavg,
                             t_wt_15_parox_runavg, t_wt_16_parox_runavg, t_wt_17_parox_runavg, t_wt_18_parox_runavg,
                             t_wt_19_parox_runavg]

gyration_runavg_wt_paroxetine = [gyr_wt_7_parox_runavg, gyr_wt_8_parox_runavg, gyr_wt_9_parox_runavg,
                                 gyr_wt_10_parox_runavg, gyr_wt_11_parox_runavg,
                                 gyr_wt_12_parox_runavg, gyr_wt_13_parox_runavg, gyr_wt_14_parox_runavg,
                                 gyr_wt_15_parox_runavg, gyr_wt_16_parox_runavg, gyr_wt_17_parox_runavg,
                                 gyr_wt_18_parox_runavg, gyr_wt_19_parox_runavg]

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

# PERCENTAGES OPEN/CLOSED STRUCTURES
df_wt_paroxetine = pd.DataFrame({'Rg': gyr_wt_parox})
print(len(df_wt_paroxetine.loc[df_wt_paroxetine['Rg'] > 2.1])/(len(df_wt_paroxetine['Rg']))*100)
print(len(df_wt_paroxetine.loc[df_wt_paroxetine['Rg'] <= 2.1])/(len(df_wt_paroxetine['Rg']))*100)

# PLOT KDE GYRATION
# plot_kde_gyration(gyr_all_wt, gyr_all_k, gyration_wt_paroxetine, gyration_k_paroxetine, figname='gyration_wt_k_parox')