from remd_distributions import *
# %%============================================================================
#                               Distance
# =============================================================================

# Read distance

# DISTANCE 56-98
dist_56_98_all_wt = read_distance_data("dist_56_98_wt.dat")
dist_56_98_all_k = read_distance_data("dist_56_98_k.dat")

# DISTANCE 17-98
dist_17_98_all_wt = read_distance_data("dist_17_98_wt.dat")
dist_17_98_all_k = read_distance_data("dist_17_98_k.dat")

# DISTANCE 48-98
dist_48_98_all_wt = read_distance_data("dist_48_98_wt.dat")
dist_48_98_all_k = read_distance_data("dist_48_98_k.dat")

# DISTANCE 44-98
dist_44_98_all_wt = read_distance_data("dist_44_98_wt.dat")
dist_44_98_all_k = read_distance_data("dist_44_98_k.dat")

# DISTANCE 86-98
dist_86_98_all_wt = read_distance_data("dist_86_98_wt.dat")
dist_86_98_all_k = read_distance_data("dist_86_98_k.dat")

# DISTANCE CYS99-CYS195 (smFRET)
dist_99_195_all_wt = read_distance_data("dist_cys99_cys195_all_wt.dat")
dist_99_195_it_wt = read_distance_data("dist_cys99_cys195_it_wt.dat")
dist_99_195_m_wt = read_distance_data("dist_cys99_cys195_m_wt.dat")
dist_99_195_trros_wt = read_distance_data("dist_cys99_cys195_trros_wt.dat")

# %% =============================================================================
#                               Plot distance
# =============================================================================

# plot_distance_distribution(2, [dist_56_98_all_wt, dist_56_98_all_k], [100, 100], ['wt', 'K141E'], ['blue', 'orange'], [0, 100], [0, 0.105], 'distance_56_98.png')

df_all_wt = pd.DataFrame({'Rg (nm)': dist_99_195_all_wt}).divide(10)
df_it_wt = pd.DataFrame({'Rg (nm)': dist_99_195_it_wt}).divide(10)
df_m_wt = pd.DataFrame({'Rg (nm)': dist_99_195_m_wt}).divide(10)
df_r_wt = pd.DataFrame({'Rg (nm)': dist_99_195_trros_wt}).divide(10)

df_wt = pd.concat([df_all_wt, df_it_wt, df_m_wt, df_r_wt], ignore_index=True)
df_wt['Homology Model'] = np.repeat(['All wt', 'I-TASSER wt', 'MODELLER wt', 'ROSETTA wt'],
                                    [int(len(dist_99_195_all_wt)), len(dist_99_195_it_wt), len(dist_99_195_m_wt),
                                     len(dist_99_195_trros_wt)])

sns.set_style('darkgrid')

fig, ax = plt.subplots(1, 1, figsize=(50, 20), gridspec_kw=dict(width_ratios=[4]))

sns.kdeplot(data=df_wt, x="Rg (nm)", hue="Homology Model", linewidth=3, shade=True, common_norm=True)

plt.axvline(6.68, color='red', linestyle='--', linewidth=4);
plt.axvline(4.94, color='red', linestyle='--', linewidth=4);

plt.text(4.54, 0.27, 'Peak 3', color='red', fontsize=60)
plt.text(6.28, 0.27, 'Peak 2', color='red', fontsize=60)

sns.set(font_scale=3.5)

ax.set_title("Distance CYS99-CYS195")
ax.set_xlabel('Distance (nm)')
ax.set_ylabel('Relative frequency')

fig.tight_layout()

fig.savefig('distance_99-195_kde_wt.png', dpi=600)
