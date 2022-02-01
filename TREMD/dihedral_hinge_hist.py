import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns


def read_dihedral_data(document_name):
  x, y = [], []
  with open(document_name) as f:
    for line in f:
        cols = line.split()

        if len(cols) == 2:
            x.append(float(cols[0]))
            y.append(float(cols[1]) % 360)
  return x, y

#90-93

time_it_wt, dihedral_it_wt_90_93_500ns = read_dihedral_data("dihedral_it_wt_90_93_268ns.dat")
time_m_wt, dihedral_m_wt_90_93_500ns = read_dihedral_data("dihedral_m_wt_90_93_207ns.dat")
time_r_wt, dihedral_r_wt_90_93_500ns = read_dihedral_data("dihedral_trros_wt_90_93_287ns.dat")

time_it_wt, dihedral_it_wt_90_93 = read_dihedral_data("dihedral_it_wt_90_93.dat")
time_m_wt, dihedral_m_wt_90_93 = read_dihedral_data("dihedral_m_wt_90_93.dat")
time_r_wt, dihedral_r_wt_90_93 = read_dihedral_data("dihedral_r_wt_90_93.dat")

time_it_k, dihedral_it_k_90_93_500ns = read_dihedral_data("dihedral_it_k_90_93_356ns.dat")
time_m_k, dihedral_m_k_90_93_500ns = read_dihedral_data("dihedral_m_k_90_93_222ns.dat")
time_r_k, dihedral_r_k_90_93_500ns = read_dihedral_data("dihedral_trros_k_90_93_287ns.dat")

time_it_k, dihedral_it_k_90_93 = read_dihedral_data("dihedral_it_k_90_93.dat")
time_m_k, dihedral_m_k_90_93 = read_dihedral_data("dihedral_m_k_90_93.dat")
time_r_k, dihedral_r_k_90_93 = read_dihedral_data("dihedral_r_k_90_93.dat")


#91-94

time_it_wt, dihedral_it_wt_91_94_500ns = read_dihedral_data("dihedral_it_wt_91_94_268ns.dat")
time_m_wt, dihedral_m_wt_91_94_500ns = read_dihedral_data("dihedral_m_wt_91_94_207ns.dat")
time_r_wt, dihedral_r_wt_91_94_500ns = read_dihedral_data("dihedral_trros_wt_91_94_287ns.dat")

time_it_wt, dihedral_it_wt_91_94 = read_dihedral_data("dihedral_it_wt_91_94.dat")
time_m_wt, dihedral_m_wt_91_94 = read_dihedral_data("dihedral_m_wt_91_94.dat")
time_r_wt, dihedral_r_wt_91_94 = read_dihedral_data("dihedral_r_wt_91_94.dat")


time_it_k, dihedral_it_k_91_94_500ns = read_dihedral_data("dihedral_it_k_91_94_356ns.dat")
time_m_k, dihedral_m_k_91_94_500ns = read_dihedral_data("dihedral_m_k_91_94_222ns.dat")
time_r_k, dihedral_r_k_91_94_500ns = read_dihedral_data("dihedral_trros_k_91_94_287ns.dat")

time_it_k, dihedral_it_k_91_94 = read_dihedral_data("dihedral_it_k_91_94.dat")
time_m_k, dihedral_m_k_91_94 = read_dihedral_data("dihedral_m_k_91_94.dat")
time_r_k, dihedral_r_k_91_94 = read_dihedral_data("dihedral_r_k_91_94.dat")

#92-95

time_it_wt, dihedral_it_wt_92_95_500ns = read_dihedral_data("dihedral_it_wt_92_95_268ns.dat")
time_m_wt, dihedral_m_wt_92_95_500ns = read_dihedral_data("dihedral_m_wt_92_95_207ns.dat")
time_r_wt, dihedral_r_wt_92_95_500ns = read_dihedral_data("dihedral_trros_wt_92_95_287ns.dat")

time_it_k, dihedral_it_k_92_95_500ns = read_dihedral_data("dihedral_it_k_92_95_356ns.dat")
time_m_k, dihedral_m_k_92_95_500ns = read_dihedral_data("dihedral_m_k_92_95_222ns.dat")
time_r_k, dihedral_r_k_92_95_500ns = read_dihedral_data("dihedral_trros_k_92_95_287ns.dat")


#93-96

time_it_wt, dihedral_it_wt_93_96_500ns = read_dihedral_data("dihedral_it_wt_93_96_268ns.dat")
time_m_wt, dihedral_m_wt_93_96_500ns = read_dihedral_data("dihedral_m_wt_93_96_207ns.dat")
time_r_wt, dihedral_r_wt_93_96_500ns = read_dihedral_data("dihedral_trros_wt_93_96_287ns.dat")

time_it_wt, dihedral_it_wt_93_96 = read_dihedral_data("dihedral_it_wt_93_96.dat")
time_m_wt, dihedral_m_wt_93_96 = read_dihedral_data("dihedral_m_wt_93_96.dat")
time_r_wt, dihedral_r_wt_93_96 = read_dihedral_data("dihedral_r_wt_93_96.dat")

time_it_k, dihedral_it_k_93_96_500ns = read_dihedral_data("dihedral_it_k_93_96_356ns.dat")
time_m_k, dihedral_m_k_93_96_500ns = read_dihedral_data("dihedral_m_k_93_96_222ns.dat")
time_r_k, dihedral_r_k_93_96_500ns = read_dihedral_data("dihedral_trros_k_93_96_287ns.dat")

time_it_k, dihedral_it_k_93_96 = read_dihedral_data("dihedral_it_k_93_96.dat")
time_m_k, dihedral_m_k_93_96 = read_dihedral_data("dihedral_m_k_93_96.dat")
time_r_k, dihedral_r_k_93_96 = read_dihedral_data("dihedral_r_k_93_96.dat")


def plot_distance_distribution(num, y, bins, labels, colors, xlim, ylim, figure_name):
  
  sns.set_style('darkgrid')

  fig = plt.figure()
  ax = fig.add_subplot(111)

  ax.set_title("Hinge dihedral 92-95 (wt)", fontsize=30)    
  ax.set_xlabel('Angle (degree)', fontsize=30)
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

  _ = plt.xticks(np.arange(xlim[0], xlim[1], 30))

  leg = ax.legend()

  plt.show()

  fig.savefig(figure_name)

#PLOT
dihedral_comparison = [dihedral_it_wt_92_95_500ns, dihedral_m_wt_92_95_500ns, dihedral_r_wt_92_95_500ns]
labels_dihedral_comparison = ['I-TASSER', 'MODELLER', 'ROSETTA']
colors = ['blue', 'orange', 'green']

plot_distance_distribution(3, dihedral_comparison, [100, 100, 100], labels_dihedral_comparison, colors, [0, 360], [0, 0.075], 'dihedral_histogram_92_95_wt.png')

