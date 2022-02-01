import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns


def read_sasa_data(document_name):
  x, y, z = [], [], []
  with open(document_name) as f:
    for line in f:
        cols = line.split()

        if len(cols) == 3:
            x.append(float(cols[0]))
            y.append(float(cols[1]))
            z.append(float(cols[2]))
  return x, y, z

#t , sasa_it_wt_500ns, hsasa_it_wt_500ns = read_sasa_data("sasa_tot_it_wt_remd_500ns.xvg")

#t , sasa_m_wt_500ns, hsasa_m_wt_500ns = read_sasa_data("sasa_tot_m_wt_remd_500ns.xvg")

#t , sasa_trros_wt_500ns, hsasa_trros_wt_500ns = read_sasa_data("sasa_tot_trros_wt_remd_500ns.xvg")

#t , sasa_it_k_500ns, hsasa_it_k_500ns = read_sasa_data("sasa_tot_it_k_remd_500ns.xvg")

#t , sasa_m_k_500ns, hsasa_m_k_500ns = read_sasa_data("sasa_tot_m_k_remd_500ns.xvg")

#t , sasa_trros_k_500ns, hsasa_trros_k_500ns = read_sasa_data("sasa_tot_trros_k_remd_500ns.xvg")


t , sasa_all_wt_500ns, hsasa_all_wt_500ns = read_sasa_data("sasa_tot_all_wt.xvg")

t , sasa_all_k_500ns, hsasa_all_k_500ns = read_sasa_data("sasa_tot_all_k.xvg")





def read_gyration_data(document_name):
  x, y = [], []
  with open(document_name) as f:
    for line in f:
        cols = line.split()

        if len(cols) == 5:
            x.append(float(cols[0]))
            y.append(float(cols[1]))
  return x, y

#t , gyr_it_wt_500ns = read_gyration_data("gyrate_it_wt_remd_500ns.xvg")

#t , gyr_m_wt_500ns = read_gyration_data("gyrate_m_wt_remd_500ns.xvg")

#t , gyr_trros_wt_500ns = read_gyration_data("gyrate_trros_wt_remd_500ns.xvg")

#t , gyr_it_k_500ns = read_gyration_data("gyrate_it_k_remd_500ns.xvg")

#t , gyr_m_k_500ns = read_gyration_data("gyrate_m_k_remd_500ns.xvg")

#t , gyr_trros_k_500ns = read_gyration_data("gyrate_trros_k_remd_500ns.xvg")


#t, gyr_all_wt_500ns = read_gyration_data("gyrate_all_wt.xvg")

#t, gyr_all_k_500ns = read_gyration_data("gyrate_all_k.xvg")



def calculate_sasa_frac(gyration, sasa):
  gyration_array = np.array(gyration)
  sasa_array = np.array(sasa)

  sasa_norm = gyration_array**2 *4 * np.pi

  sasa_frac = sasa_array / sasa_norm

  return sasa_frac

#sasa_frac_it_wt_500ns = calculate_sasa_frac(gyr_it_wt_500ns, sasa_it_wt_500ns)

#sasa_frac_m_wt_500ns = calculate_sasa_frac(gyr_m_wt_500ns, sasa_m_wt_500ns)

#sasa_frac_trros_wt_500ns = calculate_sasa_frac(gyr_trros_wt_500ns, sasa_trros_wt_500ns)

#sasa_frac_it_k_500ns = calculate_sasa_frac(gyr_it_k_500ns, sasa_it_k_500ns)

#sasa_frac_m_k_500ns = calculate_sasa_frac(gyr_m_k_500ns, sasa_m_k_500ns)

#sasa_frac_trros_k_500ns = calculate_sasa_frac(gyr_trros_k_500ns, sasa_trros_k_500ns)


#sasa_frac_all_wt_500ns = calculate_sasa_frac(gyr_all_wt_500ns, sasa_all_wt_500ns)

#sasa_frac_all_k_500ns = calculate_sasa_frac(gyr_all_k_500ns, sasa_all_k_500ns)



def calculate_hsasa_fraction(hsasa, sasa):
    
    '''Calculate fraction of hSASA/SASA '''

    hsasa_arr = np.array(hsasa)
    sasa_arr  = np.array(sasa)

    hsasa_frac = hsasa_arr  / sasa_arr

    return hsasa_frac

#CALCULATE hSASA fractions

#hsasa_frac_it_wt_500ns = calculate_hsasa_fraction(hsasa_it_wt_500ns, sasa_it_wt_500ns)

#hsasa_frac_m_wt_500ns = calculate_hsasa_fraction(hsasa_m_wt_500ns, sasa_m_wt_500ns)

#hsasa_frac_trros_wt_500ns = calculate_hsasa_fraction(hsasa_trros_wt_500ns, sasa_trros_wt_500ns)

#hsasa_frac_it_k_500ns = calculate_hsasa_fraction(hsasa_it_k_500ns, sasa_it_k_500ns)

#hsasa_frac_m_k_500ns = calculate_hsasa_fraction(hsasa_m_k_500ns, sasa_m_k_500ns)

#hsasa_frac_trros_k_500ns = calculate_hsasa_fraction(hsasa_trros_k_500ns, sasa_trros_k_500ns)


hsasa_frac_all_wt_500ns = calculate_hsasa_fraction(hsasa_all_wt_500ns, sasa_all_wt_500ns)

hsasa_frac_all_k_500ns = calculate_hsasa_fraction(hsasa_all_k_500ns, sasa_all_k_500ns)






def plot_sasa_distribution(num, y, bins, labels, colors, xlim, ylim, figure_name):
  
  sns.set_style('darkgrid')

  fig = plt.figure()
  ax = fig.add_subplot(111)

  ax.set_title("hSASA/SASA", fontsize=30)    
  ax.set_xlabel('hSASA/SASA', fontsize=30)
  ax.set_ylabel('Relative frequency', fontsize=30)

  plt.xticks(fontsize=30)
  plt.yticks(fontsize=30)
  plt.rc('legend',fontsize=30)

  plt.xlim(xlim[0], xlim[1])
  plt.ylim(ylim[0], ylim[1])

  for i in range(0,num):
    (counts_hist, hbins, patches) = ax.hist(y[i], bins=bins[i], density=True, stacked=True, color=colors[i], alpha=0.6, label=labels[i])
    for item in patches:
      item.set_height(item.get_height()/sum(counts_hist))

  _ = plt.xticks(np.arange(xlim[0], xlim[1], 0.02))

  leg = ax.legend()

  plt.show()

  fig.savefig(figure_name)



#PLOT
#sasa_comparison = [sasa_it_wt_500ns, sasa_m_wt_500ns, sasa_trros_wt_500ns, sasa_it_k_500ns, sasa_m_k_500ns, sasa_trros_k_500ns]
labels_sasa_comparison = ['I-TASSER wt', 'MODELLER wt', 'TRROSETTA wt', 'I-TASSER K141E', 'MODELLER K141E', 'TRROSETTA K141E']
colors = ['blue', 'orange', 'green', 'magenta', 'yellow', 'cyan']

#plot_sasa_distribution(6, sasa_comparison, [100, 100, 100, 100, 100, 100], labels_sasa_comparison, colors, [1.7, 3.1], [0, 0.075],'sasa_histogram_comparison_tremd_500ns_wt_k141e')

#plot_sasa_distribution(2, [sasa_all_wt_500ns, sasa_all_k_500ns], [100, 100], ['wt', 'K141E'], ['blue', 'orange'] , [110, 170], [0, 0.05],'sasa_histogram_comparison_tremd_500ns_all_wt_k141e')

#plot_sasa_distribution(2, [hsasa_all_wt_500ns, hsasa_all_k_500ns], [100, 100], ['wt', 'K141E'], ['blue', 'orange'] , [30, 65], [0, 0.05],'hsasa_histogram_all_wt_k141e')

plot_sasa_distribution(2, [hsasa_frac_all_wt_500ns, hsasa_frac_all_k_500ns], [100, 100], ['wt', 'K141E'], ['blue', 'orange'] , [0.25, 0.49], [0, 0.05],'hsasa_frac_histogram_all_wt_k141e')


def cumsum_plot(data, title, xlabel, vline, figure_name):

    #calculate cumsum
    fig= plt.figure()
    ax = fig.add_subplot(111)

    (counts_hist_1, bins_1, patches_1) = ax.hist(data[0], bins=100, density=True, stacked=True)
    density_1 = np.array(counts_hist_1 / (sum(counts_hist_1) * np.diff(bins_1)))

    cumsum_1 = np.cumsum(np.array(density_1 * np.diff(bins_1)))

    (counts_hist_2, bins_2, patches_2) = ax.hist(data[1], bins=100, density=True, stacked=True)
    density_2 = np.array(counts_hist_2 / (sum(counts_hist_2) * np.diff(bins_2)))

    cumsum_2 = np.cumsum(np.array(density_2 * np.diff(bins_2)))

    plt.close()

    #plot cumsum
    sns.set_style('darkgrid')

    fig = plt.figure()
    ax = fig.add_subplot(111)

    ax.set_title(title, fontsize=30)
    ax.set_xlabel(xlabel, fontsize=30)
    ax.set_ylabel('Cumulative frequency', fontsize=30)

    plt.xticks(fontsize=30)
    plt.yticks(fontsize=30)
    plt.rc('legend',fontsize=30)

    ax.plot(bins_1[1:], cumsum_1, linewidth=3.0, label='wt')
    ax.plot(bins_2[1:], cumsum_2, linewidth=3.0, label='K141E')

    _ = plt.xticks(np.linspace(np.minimum(bins_1[1], bins_2[1]), np.maximum(bins_1[-1], bins_2[-1]), num=8))

    if vline:
        plt.axvline(2.1, color='red', linestyle='--')

    leg = ax.legend()

    plt.show()

    fig.savefig(figure_name)


#cumsum_plot([hsasa_all_wt_500ns, hsasa_all_k_500ns], "hSASA Cumulative frequency", "hSASA ($nm^2$)", False, "hsasa_cumsum_wt_k141e")

#cumsum_plot([sasa_all_wt_500ns, sasa_all_k_500ns], "SASA Cumulative frequency", "SASA ($nm^2$)", False, "sasa_cumsum_wt_k141e")

#cumsum_plot([hsasa_frac_all_wt_500ns, hsasa_frac_all_k_500ns], "hSASA/SASA Cumulative frequency", "hSASA/SASA", False, "hsasa_frac_cumsum_wt_k141e")



