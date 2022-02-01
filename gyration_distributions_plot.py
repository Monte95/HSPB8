import numpy as np
import pandas as pd


def read_gyration_data(document_name):
    x, y = [], []
    with open(document_name) as f:
        for line in f:
            cols = line.split()

            if len(cols) == 5:
                x.append(float(cols[0]))
                y.append(float(cols[1]))
    return np.array(x), np.array(y)


def plot_kde_cumulative_seaborn(data1, data2, figname, name=''):
    df_wt = pd.DataFrame({'Rg': data1})
    df_k = pd.DataFrame({'Rg': data2})

    df_wt_k = pd.concat([df_wt, df_k], ignore_index=True)
    df_wt_k['Variant'] = np.repeat(['{}wt'.format(name), '{}K141E'.format(name)], int(len(data1)))

    sns.set_style('darkgrid')

    fig, ax = plt.subplots(1, figsize=(50, 20), gridspec_kw=dict(width_ratios=[4]))

    sns.kdeplot(data=df_wt_k, x="Rg", hue="Variant", bw_adjust=0.8, linewidth=3, shade=True, common_norm=True,
                common_grid=True, ax=ax)

    sns.set(font_scale=3)

    ax.set_title("Gyration radius", fontsize=30)
    ax.set_xlabel('Gyration radius (nm)', fontsize=30)
    ax.set_ylabel('Relative frequency', fontsize=30)

    plt.show()

    fig.savefig('{}.png'.format(figname), dpi=320)


t, gyr_all_wt_500ns = read_gyration_data("gyrate_all_wt.xvg")
t, gyr_all_k_500ns = read_gyration_data("gyrate_all_k.xvg")

plot_kde_cumulative_seaborn(gyr_all_wt_500ns, gyr_all_k_500ns, figname='gyrate_kde_wt_k', name='')
