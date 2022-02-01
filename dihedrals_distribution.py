from remd_distributions import *
# %%============================================================================
#                                   Dihedrals
# =============================================================================

# Read dihedral data

time_all_wt, dihedral_all_wt_90_93 = read_dihedral_data("dihedral_hinge_all_wt_90_93.dat")
time_all_k, dihedral_all_k_90_93 = read_dihedral_data("dihedral_hinge_all_k_90_93.dat")

time_all_wt, dihedral_all_wt_91_94 = read_dihedral_data("dihedral_hinge_all_wt_91_94.dat")
time_all_k, dihedral_all_k_91_94 = read_dihedral_data("dihedral_hinge_all_k_91_94.dat")

time_all_wt, dihedral_all_wt_92_95 = read_dihedral_data("dihedral_hinge_all_wt_92_95.dat")
time_all_k, dihedral_all_k_92_95 = read_dihedral_data("dihedral_hinge_all_k_92_95.dat")

time_all_wt, dihedral_all_wt_93_96 = read_dihedral_data("dihedral_hinge_all_wt_93_96.dat")
time_all_k, dihedral_all_k_93_96 = read_dihedral_data("dihedral_hinge_all_k_93_96.dat")

# Read gyration data
t , gyr_all_wt_500ns = read_gyration_data("gyrate_all_wt.xvg")

t , gyr_all_k_500ns = read_gyration_data("gyrate_all_k.xvg")

# %%============================================================================
#                               Plot dihedrals
# =============================================================================

# plot_dihedral_distribution(2, [dihedral_all_wt_93_96, dihedral_all_k_93_96], [100, 100],
#                            ['wt', 'K141E'], ['blue', 'orange'], [0, 190], [0, 0.075],
#                            'dihedral_hinge_all_wt_93_96.png')

# %%============================================================================
#                       Dihedrals KDE+cumulative
# =============================================================================

# plot_distribution_cumulative_seaborn(dihedral_all_wt_91_94, dihedral_all_k_91_94, gyr_all_wt_500ns, gyr_all_k_500ns,
#                                      figname='dihedrals_91_94_kde_all_wt_k',
#                                      dihedrals_name=' 91-94', name='')

# %%============================================================================
#                       Retrieve KDE data
# =============================================================================

df_wt = pd.DataFrame({'dihedrals': dihedral_all_wt_91_94[:-1], 'Rg': gyr_all_wt_500ns})
df_k = pd.DataFrame({'dihedrals': dihedral_all_k_91_94[:-1], 'Rg': gyr_all_k_500ns})

df_wt_k = pd.concat([df_wt, df_k], ignore_index=True)
df_wt_k['Variant'] = np.repeat(['wt', 'K141E'], int(len(dihedral_all_wt_91_94[:-1])))

df_wt_k.loc[df_wt_k['Rg'] <= 2.1, 'Conformation'] = 'Closed'
df_wt_k.loc[df_wt_k['Rg'] > 2.1, 'Conformation'] = 'Open'

sns.set_style('darkgrid')

fig, ax = plt.subplots(1, 1, figsize=(50, 20), gridspec_kw=dict(width_ratios=[4]))

ax = sns.kdeplot(data=df_wt_k, x="dihedrals", hue="Conformation", clip=(0.0, 180.0), bw_adjust=0.8)

plt.show()

x = ax.lines[0].get_xdata()  # Get the x data of the distribution
y = ax.lines[0].get_ydata()  # Get the y data of the distribution
peaks, _ = np.array(find_peaks(y))

for p in range(len(peaks)):
    if y[peaks[p]] >= 0.001:
        print(x[peaks[p]])

# %%============================================================================
#                           Scatterplot 2D/KDEplot dihedrals
# =============================================================================

# print(scipy.stats.pearsonr(gyr_all_wt_500ns, dist_48_98_all_wt[:-1]))
# print(scipy.stats.pearsonr(gyr_all_k_500ns, dist_48_98_all_k[:-1]))

# df_wt = pd.DataFrame({'dihedrals': dihedral_all_wt_91_94[:-1], 'Rg': gyr_all_wt_500ns})
# df_k = pd.DataFrame({'dihedrals': dihedral_all_k_91_94[:-1], 'Rg': gyr_all_k_500ns})
#
# df_wt_k = pd.concat([df_wt, df_k], ignore_index=True)
# df_wt_k['Variant'] = np.repeat(['wt', 'K141E'], int(len(dihedral_all_wt_91_94[:-1])))
#
# df_wt_k.loc[df_wt_k['Rg'] <= 2.1, 'Conformation'] = 'Closed'
# df_wt_k.loc[df_wt_k['Rg'] > 2.1, 'Conformation'] = 'Open'
#
# df_wt_k['ConformationInt'] = (df_wt_k['Conformation'] == 'Open').astype(int)
#
# print(df_wt_k.head())
# print(df_wt_k.corr())

# sns.set_style('darkgrid')
#
# sns.set(font_scale=2)
#
# fig, ax = plt.subplots(1, 1, figsize=(50, 20), gridspec_kw=dict(width_ratios=[4]))
#
# a = sns.scatterplot(data=df_wt_k, x="Rg", y='dihedrals', hue="Conformation")
# ax = sns.kdeplot(data=df_wt_k, x="Rg", y='dihedrals', hue="Conformation", levels=50, clip=(0.0, 180.0), bw_adjust=0.8)
#
# ax.set_title("Scatterplot Rg-Dihedrals 91-94", fontsize=30)
# ax.set_xlabel('Rg (nm)', fontsize=30)
# ax.set_ylabel('Dihedral angle 91-94 ($\mathrm{\AA}$)', fontsize=30)
#
# plt.show()
#
# x = ax.lines[0].get_xdata() # Get the x data of the distribution
# y = ax.lines[0].get_ydata() # Get the y data of the distribution
# peaks, _ = np.array(find_peaks(y))
#
# fig.savefig("Scatter_Rg_dihedral_91-94.png", dpi=320)
#
#
# %%============================================================================
#                               Scatterplot 3D dihedrals
# =============================================================================
#
# df_closed = df_wt_k.loc[df_wt_k['Rg'] <= 2.1, ['dihedrals', 'Rg']]
# df_open = df_wt_k.loc[df_wt_k['Rg'] > 2.1, ['dihedrals', 'Rg']]
#
# estimate kernel density of data
# kde_open = sps.gaussian_kde(df_open.values.T)
# kde_closed = sps.gaussian_kde(df_closed.values.T)
#
# get a regular grid of points over our region of interest
# xx_open, yy_open = np.meshgrid(
#     np.linspace(0, 180, 500),
#     np.linspace(2.1, 2.9, 500))
#
# xx_closed, yy_closed = np.meshgrid(
#     np.linspace(0, 180, 500),
#     np.linspace(1.65, 2.1, 500))
#
# calculate probability density on these points
# z_open = kde_open.pdf([xx_open.ravel(), yy_open.ravel()]).reshape(xx_open.shape)
# z_closed = kde_closed.pdf([xx_closed.ravel(), yy_closed.ravel()]).reshape(xx_closed.shape)
#
# plt.figure(figsize=[20, 10])
# ax = plt.axes(projection='3d')

# ax.azim = 0
# ax.elev = 40
#
# ax.scatter3D(xx_open, yy_open, z_open, c=z_open, cmap='Blues');
# ax.scatter3D(xx_closed, yy_closed, z_closed, c=z_closed, cmap='Oranges');
#
# ax.set_title('Dihedral 91-94 Open/Closed conformations')
# ax.set_xlabel('Dihedral angle')
# ax.set_ylabel('Radius of gyration (nm)')
# ax.set_zlabel('Relative frequency')
#
# %%============================================================================
#                             Scatterplot dihedrals-Rg
# =============================================================================

# dict_wt = {'Rg': gyr_all_wt_500ns,
#            'Dihedral_90_93': dihedral_all_wt_90_93[:-1],
#            'Dihedral_91_94': dihedral_all_wt_91_94[:-1],
#            'Dihedral_92_95': dihedral_all_wt_92_95[:-1],
#            'Dihedral_93_96': dihedral_all_wt_93_96[:-1]}
#
# dict_k = {'Rg': gyr_all_k_500ns,
#           'Dihedral_90_93': dihedral_all_k_90_93[:-1],
#           'Dihedral_91_94': dihedral_all_k_91_94[:-1],
#           'Dihedral_92_95': dihedral_all_k_92_95[:-1],
#           'Dihedral_93_96': dihedral_all_k_93_96[:-1]}
#
# df_wt = pd.DataFrame(data=dict_wt)
# df_k = pd.DataFrame(data=dict_k)
#
# df_wt_k = pd.concat([df_wt, df_k], ignore_index=True)
# df_wt_k['Type'] = np.repeat(['wt', 'K141E'], int(len(gyr_all_wt_500ns)))
#
# print(df_wt.corr())
# print(df_k.corr())
# print(df_wt_k.corr())
#
# fig, ax = plt.subplots()
# ax.set_title("Dihedral hinge 93-96/Rg")

# sns.scatterplot(data=df_wt_k, x='Rg', y='Dihedral_90_93', hue='Type', ax=ax)
# sns.scatterplot(data=df_wt_k, x='Rg', y='Dihedral_91_94', hue='Type')
# sns.scatterplot(data=df_wt_k, x='Rg', y='Dihedral_92_95', hue='Type')
# sns.scatterplot(data=df_wt_k, x='Rg', y='Dihedral_93_96', hue='Type')

# plt.show()

# fig.savefig("Scatter_dihedral_hinge_93_96.png", dpi=320)
