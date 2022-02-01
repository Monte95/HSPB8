import numpy as np
from remd_distributions import *

t, sasa_all_wt_500ns, hsasa_all_wt_500ns = read_sasa_data("sasa_tot_all_wt.xvg")

t, sasa_all_k_500ns, hsasa_all_k_500ns = read_sasa_data("sasa_tot_all_k.xvg")

hsasa_frac_all_wt_500ns = calculate_hsasa_fraction(hsasa_all_wt_500ns, sasa_all_wt_500ns)

hsasa_frac_all_k_500ns = calculate_hsasa_fraction(hsasa_all_k_500ns, sasa_all_k_500ns)

# plot_kde_cumulative_seaborn(hsasa_all_wt_500ns, hsasa_all_k_500ns, figname='hsasa_kde_all_wt_k',
#                             type='SASA', name='')

df_wt = pd.DataFrame({'SASA': sasa_all_wt_500ns, 'hSASA': hsasa_all_wt_500ns})
df_k = pd.DataFrame({'SASA': sasa_all_k_500ns, 'hSASA': hsasa_all_k_500ns})

df_wt_k = pd.concat([df_wt, df_k], ignore_index=True)
df_wt_k['Variant'] = np.repeat(['wt', 'K141E'], int(len(sasa_all_wt_500ns)))

sns.set_style('darkgrid')

fig, ax = plt.subplots(1, 2, figsize=(50, 20), gridspec_kw=dict(width_ratios=[4, 4]))

a = sns.kdeplot(data=df_wt_k, x="SASA", hue="Variant", bw_adjust=0.8, linewidth=3, shade=True, common_norm=True,
                common_grid=True, ax=ax[0])

ax[0].set_title("SASA", fontsize=20)
ax[0].set_xlabel('SASA ($\mathrm{nm^2}$)', fontsize=20)
ax[0].set_ylabel('Relative frequency', fontsize=20)

b = sns.kdeplot(data=df_wt_k, x="hSASA", hue="Variant", bw_adjust=0.8, linewidth=3, shade=True, common_norm=True,
                common_grid=True, ax=ax[1])

print(np.median(hsasa_all_k_500ns))
# plt.axvline([np.median(sasa_all_wt_500ns)], color='blue', linestyle='--', linewidth=2)
# plt.axvline([np.median(sasa_all_k_500ns)], color='orange', linestyle='--', linewidth=2)

ax[1].set_title("hSASA", fontsize=20)
ax[1].set_xlabel('hSASA ($\mathrm{nm^2}$)', fontsize=20)
ax[1].set_ylabel('Relative frequency', fontsize=20)

for label in (ax[0].get_xticklabels() + ax[0].get_yticklabels()):
    label.set_fontsize(12)

for label in (ax[1].get_xticklabels() + ax[1].get_yticklabels()):
    label.set_fontsize(12)

# fig.savefig('sasa_hsasa_kde_all_wt_k.png', dpi=320)

# plt.show()
