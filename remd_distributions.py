""#%%============================================================================
#                           Import Libraries
# =============================================================================

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from scipy.signal import find_peaks
import scipy.stats as sps

#%%============================================================================
#                              Functions
# =============================================================================


def read_gyration_data(document_name):
  x, y = [], []
  with open(document_name) as f:
    for line in f:
        cols = line.split()

        if len(cols) == 5:
            x.append(float(cols[0]))
            y.append(float(cols[1]))
  return np.array(x), np.array(y)


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


def read_dihedral_data(document_name):
  x, y = [], []
  with open(document_name) as f:
    for line in f:
        cols = line.split()

        if len(cols) == 2:
            x.append(float(cols[0]))
            y.append(abs(float(cols[1])))
            #y.append(float(cols[1]) % 360)
            #y.append(float(cols[1]))
  return x, y


def read_rmsd_data(document_name):
  x, y = [], []
  with open(document_name) as f:
    for line in f:
        cols = line.split()

        if len(cols) == 2:
            x.append(float(cols[0]))
            y.append(float(cols[1]))
  return x, y


def read_distance_data(document_name):
  x, y = [], []
  with open(document_name) as f:
    for line in f:
        cols = line.split()

        if len(cols) == 2:
            x.append(float(cols[0]))
            y.append(float(cols[1]))
  return y


def calculate_hsasa_fraction(hsasa, sasa):
    
    '''Calculate fraction of hSASA/SASA '''

    hsasa_arr = np.array(hsasa)
    sasa_arr  = np.array(sasa)

    hsasa_frac = hsasa_arr  / sasa_arr

    return hsasa_frac


def running_avg(window_size, y, num_ns):
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

  ax.set_title("Gyration radius", fontsize=30)
  ax.set_xlabel('Time (ns)', fontsize=30)
  ax.set_ylabel('Gyration radius (nm)', fontsize=30)

  plt.xticks(fontsize=30)
  plt.yticks(fontsize=30)
  plt.rc('legend',fontsize=30)

  #plt.xlim(xlim[0], xlim[1])
  plt.ylim(ylim[0], ylim[1])

  for i in range(0,num):
    ax.plot(x[i], y[i], color=colors[i], alpha=0.6, label=labels[i], linewidth=3)

  plt.hlines(2.08, 0, 500, linestyles='dashed')


  #_ = plt.xticks(np.arange(xlim[0], xlim[1], 100000))

  leg = ax.legend()

  plt.show()

  fig.savefig(figure_name)


def plot_rmsd(num, x, y, labels, colors, ylim, figure_name):

  sns.set_style('darkgrid')

  fig = plt.figure(figsize=(50,20))
  ax = fig.add_subplot(111)

  ax.set_title("RMSD", fontsize=30)
  ax.set_xlabel('Time (ns)', fontsize=30)
  ax.set_ylabel('RMSD (nm)', fontsize=30)

  plt.xticks(fontsize=30)
  plt.yticks(fontsize=30)
  plt.rc('legend',fontsize=30)

  #plt.xlim(xlim[0], xlim[1])
  plt.ylim(ylim[0], ylim[1])

  for i in range(0,num):
    ax.plot(x[i], y[i], color=colors[i], alpha=0.6, label=labels[i], linewidth=3)

  #_ = plt.xticks(np.arange(xlim[0], xlim[1], 100000))

  ax.legend()

  plt.show()

  fig.savefig(figure_name, dpi=320)


def plot_gyration_distribution(num, y, bins, labels, colors, xlim, ylim, figure_name):

  sns.set_style('darkgrid')

  fig = plt.figure()
  ax = fig.add_subplot(111)

  ax.set_title("Gyration radius", fontsize=30)
  ax.set_xlabel('Gyration radius (nm)', fontsize=30)
  ax.set_ylabel('Relative frequency', fontsize=30)

  plt.xticks(fontsize=30)
  plt.yticks(fontsize=30)
  plt.rc('legend',fontsize=20)

  plt.xlim(xlim[0], xlim[1])
  plt.ylim(ylim[0], ylim[1])

  for i in range(0,num):
    (counts_hist, hbins, patches) = ax.hist(y[i], bins=bins[i], density=True, stacked=True, color=colors[i], alpha=0.6, label=labels[i])
    for item in patches:
      item.set_height(item.get_height()/sum(counts_hist))

  _ = plt.xticks(np.arange(xlim[0], xlim[1], 0.1))
  
  plt.axvline(2.1, color='red', linestyle='--')
  #plt.hlines(0.003, 0, 500, color='red', linestyles='--')

  leg = ax.legend()

  plt.show()

  fig.savefig(figure_name)


def plot_dihedral_distribution(num, y, bins, labels, colors, xlim, ylim, figure_name):

  sns.set_style('darkgrid')

  fig = plt.figure()
  ax = fig.add_subplot(111)

  ax.set_title("Hinge dihedral 93-96 (wt)", fontsize=30)
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


def plot_distance_distribution(num, y, bins, labels, colors, xlim, ylim, figure_name):

  sns.set_style('darkgrid')

  fig = plt.figure()
  ax = fig.add_subplot(111)

  ax.set_title("Distance residues 56-98", fontsize=30)
  ax.set_xlabel('Distance (A)', fontsize=30)
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
  

def plot_kde_cumulative_seaborn(data1, data2, figname, type='Rg', name=''):
    
    if type == 'Rg':
        
        df_wt = pd.DataFrame({'Rg' :data1})
        df_k  = pd.DataFrame({'Rg' :data2})
         
        df_wt_k = pd.concat([df_wt, df_k], ignore_index=True)
        df_wt_k['Variant'] = np.repeat(['{}wt'.format(name), '{}K141E'.format(name)], int(len(data1)))

        sns.set_style('darkgrid')

        fig, ax = plt.subplots(1,2,figsize=(50,20), gridspec_kw=dict(width_ratios=[4,4]))

        plt.axvline(2.1, color='red', linestyle='--', linewidth=4);

        a = sns.kdeplot(data=df_wt_k, x="Rg", hue="Variant", bw_adjust=0.8, linewidth=3, shade=True, common_norm=True, common_grid=True, ax=ax[0])
                
        ax[0].text(1.85, 1.3, 'Closed', color='red', fontsize=30)
        ax[0].text(2.66, 1.3, 'Open', color='red', fontsize=30)
        ax[0].set_ylim(0, 1.7)

        sns.set(font_scale = 3)

        ax[0].set_title("Gyration radius", fontsize=40)
        ax[0].set_xlabel('Gyration radius (nm)', fontsize=40)
        ax[0].set_ylabel('Relative frequency', fontsize=40)

        b = sns.ecdfplot(data=df_wt_k, x="Rg", hue="Variant", linewidth=3, ax=ax[1])
        
        ax[1].text(1.86, 0.925, 'Closed', color='red', fontsize=30)
        ax[1].text(2.67, 0.925, 'Open', color='red', fontsize=30)

        ax[1].set_title("Cumulative Gyration radius", fontsize=40)
        ax[1].set_xlabel('Gyration radius (nm)', fontsize=40)
        ax[1].set_ylabel('Cumulative frequency', fontsize=40)

        ax[0].axvline(2.1, color='red', linestyle='--', linewidth=4);

        plt.rc('legend', fontsize=30)


        fig.tight_layout()

        fig.savefig('{}.png'.format(figname), dpi=320)

    elif type == 'Rg_nocumsum':

        df_wt = pd.DataFrame({'Rg': data1})
        df_k = pd.DataFrame({'Rg': data2})

        df_wt_k = pd.concat([df_wt, df_k], ignore_index=True)
        df_wt_k['Variant'] = np.repeat(['{}wt'.format(name), '{}K141E'.format(name)], int(len(data1)))

        sns.set_style('darkgrid')

        fig, ax = plt.subplots(1, figsize=(50, 20), gridspec_kw=dict(width_ratios=[4]))

        sns.kdeplot(data=df_wt_k, x="Rg", hue="Variant", bw_adjust=0.8, linewidth=3, shade=True, common_norm=True,
                    common_grid=True, ax=ax)

        sns.set(font_scale=3)

        # plt.axvline([2.04], color='blue', linestyle='--', linewidth=2)
        # plt.axvline([2.41], color='blue', linestyle='--', linewidth=2)
        # plt.axvline([2.30], color='blue', linestyle='--', linewidth=2)
        # plt.axvline([2.36], color='blue', linestyle='--', linewidth=2)
        # plt.axvline([2.31], color='blue', linestyle='--', linewidth=2)
        # plt.axvline([2.46], color='blue', linestyle='--', linewidth=2)
        # plt.axvline([2.21], color='blue', linestyle='--', linewidth=2)
        # plt.axvline([2.25], color='blue', linestyle='--', linewidth=2)
        # plt.axvline([1.90], color='blue', linestyle='--', linewidth=2)
        # plt.axvline([2.09], color='blue', linestyle='--', linewidth=2)
        # plt.axvline([2.71], color='blue', linestyle='--', linewidth=2)
        # plt.axvline([2.27], color='blue', linestyle='--', linewidth=2)
        # plt.axvline([2.32], color='blue', linestyle='--', linewidth=2)
        # plt.axvline([1.94], color='orange', linestyle='--', linewidth=2)
        # plt.axvline([2.00], color='orange', linestyle='--', linewidth=2)
        # plt.axvline([1.83], color='orange', linestyle='--', linewidth=2)
        # plt.axvline([2.14], color='orange', linestyle='--', linewidth=2)
        # plt.axvline([2.67], color='orange', linestyle='--', linewidth=2)
        # plt.axvline([2.71], color='orange', linestyle='--', linewidth=2)
        # plt.axvline([2.33], color='orange', linestyle='--', linewidth=2)
        # plt.axvline([2.46], color='orange', linestyle='--', linewidth=2)
        # plt.axvline([2.04], color='orange', linestyle='--', linewidth=2)

        ax.set_title("Gyration radius", fontsize=30)
        ax.set_xlabel('Gyration radius (nm)', fontsize=30)
        ax.set_ylabel('Relative frequency', fontsize=30)

        #fig.tight_layout()

        plt.show()

        fig.savefig('{}.png'.format(figname), dpi=320)

    elif type == 'SASA':
    
        df_wt = pd.DataFrame({'SASA':data1})
        df_k = pd.DataFrame({'SASA':data2})
         
        df_wt_k = pd.concat([df_wt, df_k], ignore_index=True)
        df_wt_k['Variant'] = np.repeat(['{}wt'.format(name), '{}K141E'.format(name)], int(len(data1)))

        sns.set_style('darkgrid')

        fig, ax = plt.subplots(1,2,figsize=(50,20), gridspec_kw=dict(width_ratios=[4, 4]))

        a = sns.kdeplot(data=df_wt_k, x="SASA", hue="Variant", bw_adjust=0.8, linewidth=3, shade=True, common_norm=True,
                        common_grid=True, ax=ax[0])

        # sns.set(font_scale = 3)

        ax[0].set_title("SASA", fontsize=20)
        ax[0].set_xlabel('hSASA ($\mathrm{nm^2}$)', fontsize=20)
        ax[0].set_ylabel('Relative frequency', fontsize=20)

        b = sns.ecdfplot(data=df_wt_k, x="SASA", hue="Variant", linewidth=2, ax=ax[1])

        ax[1].set_title("Cumulative hSASA", fontsize=20)
        ax[1].set_xlabel('hSASA ($\mathrm{nm^2}$)', fontsize=20)
        ax[1].set_ylabel('Cumulative frequency', fontsize=20)
                
        # fig.tight_layout()

        fig.savefig('{}.png'.format(figname), dpi=320)
        
    elif type == 'SASA_nocumsum':
        
        df_wt = pd.DataFrame({'SASA' :data1})
        df_k  = pd.DataFrame({'SASA' :data2})
         
        df_wt_k = pd.concat([df_wt, df_k], ignore_index=True)
        df_wt_k['Variant'] = np.repeat(['{}wt'.format(name), '{}K141E'.format(name)], int(len(data1)))

        sns.set_style('darkgrid')

        fig, ax = plt.subplots(1, figsize=(50,20), gridspec_kw=dict(width_ratios=[4]))

        sns.kdeplot(data=df_wt_k, x="SASA", hue="Variant", bw_adjust=0.8, linewidth=3, shade=True, common_norm=True, common_grid=True, ax=ax)

        plt.axvline([0.332], color='blue', linestyle='--', linewidth=2)
        plt.axvline([0.333], color='blue', linestyle='--', linewidth=2)
        plt.axvline([0.327], color='blue', linestyle='--', linewidth=2)
        plt.axvline([0.343], color='blue', linestyle='--', linewidth=2)
        plt.axvline([0.373], color='blue', linestyle='--', linewidth=2)
        plt.axvline([0.346], color='blue', linestyle='--', linewidth=2)
        plt.axvline([0.386], color='blue', linestyle='--', linewidth=2)
        plt.axvline([0.384], color='blue', linestyle='--', linewidth=2)
        plt.axvline([0.317], color='blue', linestyle='--', linewidth=2)
        plt.axvline([0.316], color='blue', linestyle='--', linewidth=2)
        plt.axvline([0.341], color='blue', linestyle='--', linewidth=2)
        plt.axvline([0.316], color='blue', linestyle='--', linewidth=2)
        plt.axvline([0.334], color='blue', linestyle='--', linewidth=2)
        plt.axvline([0.370], color='orange', linestyle='--', linewidth=2)
        plt.axvline([0.383], color='orange', linestyle='--', linewidth=2)
        plt.axvline([0.386], color='orange', linestyle='--', linewidth=2)
        plt.axvline([0.310], color='orange', linestyle='--', linewidth=2)
        plt.axvline([0.347], color='orange', linestyle='--', linewidth=2)
        plt.axvline([0.354], color='orange', linestyle='--', linewidth=2)
        plt.axvline([0.345], color='orange', linestyle='--', linewidth=2)
        plt.axvline([0.338], color='orange', linestyle='--', linewidth=2)
        plt.axvline([0.326], color='orange', linestyle='--', linewidth=2)

        sns.set(font_scale=3)

        ax.set_title("hSASA/SASA", fontsize=30)
        ax.set_xlabel('hSASA/SASA', fontsize=30)
        ax.set_ylabel('Relative frequency', fontsize=30)

        #fig.tight_layout()

        plt.show()

        fig.savefig('{}.png'.format(figname), dpi=320)        


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
    
    print(bins_2)
    print(cumsum_2-cumsum_1)

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
    
    _ = plt.xticks(np.linspace(np.minimum(bins_1[1], bins_2[1]), np.maximum(bins_1[-1], bins_2[-1]), num=10))

    if vline:
        plt.axvline(2.1, color='red', linestyle='--')

    leg = ax.legend()

    plt.show()

    fig.savefig(figure_name)


def plot_distribution_cumulative_seaborn(data1, data2, Rg1, Rg2, figname, dihedrals_name='', name=''):
        
    df_wt= pd.DataFrame({'dihedrals': data1[:-1], 'Rg': Rg1})
    df_k= pd.DataFrame({'dihedrals': data2[:-1], 'Rg': Rg2})
         
    df_wt_k = pd.concat([df_wt, df_k], ignore_index=True)
    df_wt_k['Variant'] = np.repeat(['{}wt'.format(name), '{}K141E'.format(name)], int(len(data1[:-1])))
    
    df_wt_k.loc[df_wt_k['Rg'] <= 2.1 ,'Conformation'] = 'Closed'
    df_wt_k.loc[df_wt_k['Rg']  > 2.1 ,'Conformation'] = 'Open'
    
    sns.set_style('darkgrid')

    fig, ax = plt.subplots(1,2,figsize=(50,20), gridspec_kw=dict(width_ratios=[4,4]))

    a = sns.kdeplot(data=df_wt_k, x="dihedrals", hue="Conformation", clip=(0.0, 180.0), bw_adjust=0.8, linewidth=3, shade=True, common_norm=True, common_grid=True, ax=ax[0])

    sns.set(font_scale = 2)
    
    ax[0].set_title("Dihedrals{}".format(dihedrals_name), fontsize=30)
    ax[0].set_xlabel('Dihedral angle', fontsize=30)
    ax[0].set_ylabel('Relative frequency', fontsize=30)

    b = sns.ecdfplot(data=df_wt_k, x="dihedrals", hue="Conformation", linewidth=3, ax=ax[1])

    ax[1].set_title("Cumulative Dihedrals", fontsize=30)
    ax[1].set_xlabel('Dihedral angle', fontsize=30)
    ax[1].set_ylabel('Cumulative frequency', fontsize=30)

    fig.tight_layout()

    plt.show()

    # fig.savefig('{}.png'.format(figname), dpi=320)
    







