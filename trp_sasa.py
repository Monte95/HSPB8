import mdtraj
import numpy as np
import pandas as pd
from remd_distributions import read_sasa_data
from trp import *

SASA_TRP_SUPERFICIE = 0.013652

t, sasa_tot_wt, sasa_trp_wt = read_sasa_data("sasa_tot_trp_wt.xvg")
t, sasa_tot_k, sasa_trp_k = read_sasa_data("sasa_tot_trp_k.xvg")

trp_wt_sasa_df = pd.DataFrame({'Time': t,
                               'SASA': sasa_tot_wt,
                               'SASA TRP': sasa_trp_wt,
                               'SASA TRP/SASA': np.array(sasa_trp_wt) / np.array(sasa_tot_wt),
                               'Variant': np.repeat('wt', len(sasa_tot_wt))})

trp_k_sasa_df = pd.DataFrame({'Time': t,
                              'SASA': sasa_tot_k,
                              'SASA TRP': sasa_trp_k,
                              'SASA TRP/SASA': np.array(sasa_trp_k) / np.array(sasa_tot_k),
                              'Variant': np.repeat('K141E', len(sasa_tot_wt))})

trp_wt_k_sasa_df = pd.concat([trp_wt_sasa_df, trp_k_sasa_df], ignore_index=True)

t, sasa_tot_wt, sasa_trp_48_wt = read_sasa_data('sasa_tot_trp_48_wt.xvg')
t, sasa_tot_wt, sasa_trp_51_wt = read_sasa_data('sasa_tot_trp_51_wt.xvg')
t, sasa_tot_wt, sasa_trp_60_wt = read_sasa_data('sasa_tot_trp_60_wt.xvg')
t, sasa_tot_wt, sasa_trp_96_wt = read_sasa_data('sasa_tot_trp_96_wt.xvg')

t, sasa_tot_k, sasa_trp_48_k = read_sasa_data('sasa_tot_trp_48_k.xvg')
t, sasa_tot_k, sasa_trp_51_k = read_sasa_data('sasa_tot_trp_51_k.xvg')
t, sasa_tot_k, sasa_trp_60_k = read_sasa_data('sasa_tot_trp_60_k.xvg')
t, sasa_tot_k, sasa_trp_96_k = read_sasa_data('sasa_tot_trp_96_k.xvg')

trp_48_wt_sasa_df = pd.DataFrame({'Time': t,
                                  'SASA': sasa_tot_wt,
                                  'SASA TRP': sasa_trp_48_wt,
                                  'SASA TRP/SASA': np.array(sasa_trp_48_wt) / np.array(sasa_tot_wt),
                                  'Variant': np.repeat('wt', len(sasa_tot_wt)),
                                  'TRP Residue': np.repeat('48', len(sasa_tot_wt))})

trp_51_wt_sasa_df = pd.DataFrame({'Time': t,
                                  'SASA': sasa_tot_wt,
                                  'SASA TRP': sasa_trp_51_wt,
                                  'SASA TRP/SASA': np.array(sasa_trp_51_wt) / np.array(sasa_tot_wt),
                                  'Variant': np.repeat('wt', len(sasa_tot_wt)),
                                  'TRP Residue': np.repeat('51', len(sasa_tot_wt))})

trp_60_wt_sasa_df = pd.DataFrame({'Time': t,
                                  'SASA': sasa_tot_wt,
                                  'SASA TRP': sasa_trp_60_wt,
                                  'SASA TRP/SASA': np.array(sasa_trp_60_wt) / np.array(sasa_tot_wt),
                                  'Variant': np.repeat('wt', len(sasa_tot_wt)),
                                  'TRP Residue': np.repeat('60', len(sasa_tot_wt))})

trp_96_wt_sasa_df = pd.DataFrame({'Time': t,
                                  'SASA': sasa_tot_wt,
                                  'SASA TRP': sasa_trp_96_wt,
                                  'SASA TRP/SASA': np.array(sasa_trp_96_wt) / np.array(sasa_tot_wt),
                                  'Variant': np.repeat('wt', len(sasa_tot_wt)),
                                  'TRP Residue': np.repeat('96', len(sasa_tot_wt))})

##
trp_48_k_sasa_df = pd.DataFrame({'Time': t,
                                  'SASA': sasa_tot_k,
                                  'SASA TRP': sasa_trp_48_k,
                                  'SASA TRP/SASA': np.array(sasa_trp_48_k) / np.array(sasa_tot_k),
                                  'Variant': np.repeat('K141E', len(sasa_tot_k)),
                                  'TRP Residue': np.repeat('48', len(sasa_tot_k))})

trp_51_k_sasa_df = pd.DataFrame({'Time': t,
                                  'SASA': sasa_tot_k,
                                  'SASA TRP': sasa_trp_51_k,
                                  'SASA TRP/SASA': np.array(sasa_trp_51_k) / np.array(sasa_tot_k),
                                  'Variant': np.repeat('K141E', len(sasa_tot_k)),
                                  'TRP Residue': np.repeat('51', len(sasa_tot_k))})

trp_60_k_sasa_df = pd.DataFrame({'Time': t,
                                  'SASA': sasa_tot_k,
                                  'SASA TRP': sasa_trp_60_k,
                                  'SASA TRP/SASA': np.array(sasa_trp_60_k) / np.array(sasa_tot_k),
                                  'Variant': np.repeat('K141E', len(sasa_tot_k)),
                                  'TRP Residue': np.repeat('60', len(sasa_tot_k))})

trp_96_k_sasa_df = pd.DataFrame({'Time': t,
                                  'SASA': sasa_tot_k,
                                  'SASA TRP': sasa_trp_96_k,
                                  'SASA TRP/SASA': np.array(sasa_trp_96_k) / np.array(sasa_tot_k),
                                  'Variant': np.repeat('K141E', len(sasa_tot_k)),
                                  'TRP Residue': np.repeat('96', len(sasa_tot_k))})

trp_singoli_wt_k_sasa_df = pd.concat([trp_48_wt_sasa_df, trp_51_wt_sasa_df, trp_60_wt_sasa_df, trp_96_wt_sasa_df,
                                      trp_48_k_sasa_df, trp_51_k_sasa_df, trp_60_k_sasa_df, trp_96_k_sasa_df],
                                     ignore_index=True)

# print(trp_51_wt_sasa_df.loc[trp_51_wt_sasa_df['SASA TRP/SASA'] == max(trp_51_wt_sasa_df['SASA TRP/SASA'])])
print()

# PLOT
sns.set_style('darkgrid')

fig, ax = plt.subplots(1, figsize=(50, 20), gridspec_kw=dict(width_ratios=[4]))

# sns.violinplot(data=trp_singoli_wt_k_sasa_df, x="TRP Residue", y='SASA TRP/SASA', hue="Variant", split=True)

sns.boxplot(data=trp_singoli_wt_k_sasa_df, x='TRP Residue', y='SASA TRP/SASA', hue='Variant', linewidth=2.5)

#sns.set(font_scale=3)
ax.set_title("TRP SASA", fontsize=20)
ax.set_xlabel('TRP Residue', fontsize=20)
ax.set_ylabel('SASA TRP/SASA', fontsize=20)

for label in (ax.get_xticklabels() + ax.get_yticklabels()):
    label.set_fontsize(14)

leg = ax.legend(prop={"size": 20})


plt.show()

fig.savefig('boxplot_trp_sasa_wt_k.png', dpi=320)
