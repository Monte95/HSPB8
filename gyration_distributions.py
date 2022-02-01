from remd_distributions import *

# %%============================================================================
#                           Gyration Radius
# =============================================================================

# READ DATA

t_it_wt , gyr_it_wt_500ns = read_gyration_data("gyrate_it_wt_remd_500ns.xvg")

t_m_wt , gyr_m_wt_500ns = read_gyration_data("gyrate_m_wt_remd_500ns.xvg")

t_trros_wt , gyr_trros_wt_500ns = read_gyration_data("gyrate_trros_wt_remd_500ns.xvg")

t_it_k , gyr_it_k_500ns = read_gyration_data("gyrate_it_k_remd_500ns.xvg")

t_m_k , gyr_m_k_500ns = read_gyration_data("gyrate_m_k_remd_500ns.xvg")

t_trros_k , gyr_trros_k_500ns = read_gyration_data("gyrate_trros_k_remd_500ns.xvg")

# t_it_wt_md , gyr_it_wt_md = read_gyration_data("gyrate_it_wt_md_nter.xvg")

# t_m_wt_md , gyr_m_wt_md = read_gyration_data("gyrate_m_wt_md_nter.xvg")

# t_trros_wt_md , gyr_trros_wt_md = read_gyration_data("gyrate_trros_wt_md_nter.xvg")

# t_it_k_md , gyr_it_k_md = read_gyration_data("gyrate_it_k_md_nter.xvg")

# t_m_k_md , gyr_m_k_md = read_gyration_data("gyrate_m_k_md_nter.xvg")

# t_trros_k_md , gyr_trros_k_md = read_gyration_data("gyrate_trros_k_md_nter.xvg")


t , gyr_all_wt_500ns = read_gyration_data("gyrate_all_wt.xvg")

t , gyr_all_k_500ns = read_gyration_data("gyrate_all_k.xvg")



# %% =============================================================================
#                           Running averages
# =============================================================================

window_size = 100

t_it_wt_runavg, gyr_it_wt_500ns_runavg = running_avg(window_size, gyr_it_wt_500ns, len(t_it_wt ) /10)

t_it_k_runavg, gyr_it_k_500ns_runavg = running_avg(window_size, gyr_it_k_500ns, len(t_it_k ) /10)

t_m_wt_runavg, gyr_m_wt_500ns_runavg = running_avg(window_size, gyr_m_wt_500ns, len(t_m_wt ) /10)

t_m_k_runavg, gyr_m_k_500ns_runavg = running_avg(window_size, gyr_m_k_500ns, len(t_m_k ) /10)

t_trros_wt_runavg, gyr_trros_wt_500ns_runavg = running_avg(window_size, gyr_trros_wt_500ns, len(t_trros_wt ) /10)

t_trros_k_runavg, gyr_trros_k_500ns_runavg = running_avg(window_size, gyr_trros_k_500ns, len(t_trros_k ) /10)

# t_it_wt_md_runavg, gyr_it_wt_md_runavg = running_avg(window_size, gyr_it_wt_md, len(t_it_wt_md)/10)

# t_m_wt_md_runavg, gyr_m_wt_md_runavg = running_avg(window_size, gyr_m_wt_md, len(t_m_wt_md)/10)

# t_trros_wt_md_runavg, gyr_trros_wt_md_runavg = running_avg(window_size, gyr_trros_wt_md, len(t_trros_wt_md)/10)

# t_it_k_md_runavg, gyr_it_k_md_runavg = running_avg(window_size, gyr_it_k_md, len(t_it_k_md)/10)

# t_m_k_md_runavg, gyr_m_k_md_runavg = running_avg(window_size, gyr_m_k_md, len(t_m_k_md)/10)

# t_trros_k_md_runavg, gyr_trros_k_md_runavg = running_avg(window_size, gyr_trros_k_md, len(t_trros_k_md)/10)

t_all_wt_runavg, gyr_all_wt_runavg = running_avg(window_size, gyr_all_wt_500ns, len(gyr_all_wt_500ns) /10)
t_all_k_runavg, gyr_all_k_runavg = running_avg(window_size, gyr_all_k_500ns, len(gyr_all_k_500ns) /10)

# COUNT OPEN/CLOSED AND CLOSED/OPEN TRANSITIONS
count_oc_wt = 0
count_co_wt = 0

count_oc_k = 0
count_co_k = 0

gyration_runavg_wt = [gyr_it_wt_500ns_runavg, gyr_m_wt_500ns_runavg, gyr_trros_wt_500ns_runavg]
gyration_runavg_k = [gyr_it_k_500ns_runavg, gyr_m_k_500ns_runavg, gyr_trros_k_500ns_runavg]

gyration_wt = [gyr_it_wt_500ns, gyr_m_wt_500ns, gyr_trros_wt_500ns]
gyration_k = [gyr_it_k_500ns, gyr_m_k_500ns, gyr_trros_k_500ns]

for gyr_wt in gyration_wt:
    for j in range(1, len(gyr_wt)):
        if gyr_wt[j] < 2.1 and gyr_wt[j-1] > 2.1:
            count_oc_wt += 1

        if gyr_wt[j] > 2.1 and gyr_wt[j-1] < 2.1:
            count_co_wt += 1

for gyr_k in gyration_k:
    for j in range(1, len(gyr_k)):
        if gyr_k[j] < 2.1 and gyr_k[j-1] > 2.1:
            count_oc_k += 1

        if gyr_k[j] > 2.1 and gyr_k[j-1] < 2.1:
            count_co_k += 1

print('closed -> open wt: {}'.format(count_co_wt))
print('open -> closed wt: {}'.format(count_oc_wt))
print('closed -> open k: {}'.format(count_co_k))
print('open -> closed k: {}'.format(count_oc_k))

# %%============================================================================
#                                   Plot gyration
# =============================================================================

# time_runavg = [t_it_wt_runavg, t_m_wt_runavg, t_trros_wt_runavg, t_it_k_runavg, t_m_k_runavg, t_trros_k_runavg]
# time_md_runavg = [t_it_wt_md_runavg, t_m_wt_md_runavg, t_trros_wt_md_runavg, t_it_k_md_runavg, t_m_k_md_runavg, t_trros_k_md_runavg]
# time = [t_it_wt, t_m_wt, t_trros_wt, t_it_k, t_m_k, t_trros_k]

# gyration_runavg = [gyr_it_wt_500ns_runavg, gyr_m_wt_500ns_runavg, gyr_trros_wt_500ns_runavg, gyr_it_k_500ns_runavg, gyr_m_k_500ns_runavg, gyr_trros_k_500ns_runavg]
# gyration_md_runavg = [gyr_it_wt_md_runavg, gyr_m_wt_md_runavg, gyr_trros_wt_md_runavg, gyr_it_k_md_runavg, gyr_m_k_md_runavg, gyr_trros_k_md_runavg]
# gyration_comparison = [gyr_it_wt_500ns, gyr_m_wt_500ns, gyr_trros_wt_500ns, gyr_it_k_500ns, gyr_m_k_500ns, gyr_trros_k_500ns]

labels_gyration_comparison = ['I-TASSER wt', 'MODELLER wt', 'TRROSETTA wt', 'I-TASSER K141E', 'MODELLER K141E', 'TRROSETTA K141E']
colors = ['blue', 'orange', 'green', 'magenta', 'yellow', 'cyan']

# plot_gyration_distribution(6, gyration_comparison, [100, 100, 100, 100, 100, 100], labels_gyration_comparison, colors, [1.7, 3.1], [0, 0.075],'gyration_histogram_comparison_tremd_500ns_wt_k141e')
# plot_gyration_distribution(4, [gyr_trros_wt_500ns, gyr_m_wt_500ns, gyr_it_wt_500ns, gyr_all_wt_500ns], [100, 100, 100, 100], ['TRROSETTA wt', 'MODELLER wt', 'I-TASSER wt', 'All HM wt'], ['blue', 'orange', 'green', 'magenta'] , [1.7, 3.1], [0, 0.075],'prova_gyration')

# plot_gyration_distribution(2, [gyr_all_wt_500ns, gyr_all_k_500ns], [70, 70], ['wt', 'K141E'], ['blue', 'orange'] , [1.7, 3.1], [0, 0.05],'gyration_all_wt_k141e')

# plot_gyration(6, time_runavg, gyration_runavg, labels_gyration_comparison, colors, [0, 3.3], 'gyration_time_runavg')

# plot_gyration(6, time, gyration_comparison, labels_gyration_comparison, colors, [0, 3.2], 'gyration_time_tremd500ns')


# %%============================================================================
#                                    Read SASA
# =============================================================================

_ , sasa_it_wt_500ns, hsasa_it_wt_500ns = read_sasa_data("sasa_tot_it_wt_remd_500ns.xvg")

_ , sasa_m_wt_500ns, hsasa_m_wt_500ns = read_sasa_data("sasa_tot_m_wt_remd_500ns.xvg")

_ , sasa_trros_wt_500ns, hsasa_trros_wt_500ns = read_sasa_data("sasa_tot_trros_wt_remd_500ns.xvg")

_ , sasa_it_k_500ns, hsasa_it_k_500ns = read_sasa_data("sasa_tot_it_k_remd_500ns.xvg")

_ , sasa_m_k_500ns, hsasa_m_k_500ns = read_sasa_data("sasa_tot_m_k_remd_500ns.xvg")

_ , sasa_trros_k_500ns, hsasa_trros_k_500ns = read_sasa_data("sasa_tot_trros_k_remd_500ns.xvg")


t , sasa_all_wt_500ns, hsasa_all_wt_500ns = read_sasa_data("sasa_tot_all_wt.xvg")

t , sasa_all_k_500ns, hsasa_all_k_500ns = read_sasa_data("sasa_tot_all_k.xvg")


hsasa_frac_it_wt_500ns = calculate_hsasa_fraction(hsasa_it_wt_500ns, sasa_it_wt_500ns)

hsasa_frac_m_wt_500ns = calculate_hsasa_fraction(hsasa_m_wt_500ns, sasa_m_wt_500ns)

hsasa_frac_trros_wt_500ns = calculate_hsasa_fraction(hsasa_trros_wt_500ns, sasa_trros_wt_500ns)

hsasa_frac_it_k_500ns = calculate_hsasa_fraction(hsasa_it_k_500ns, sasa_it_k_500ns)

hsasa_frac_m_k_500ns = calculate_hsasa_fraction(hsasa_m_k_500ns, sasa_m_k_500ns)

hsasa_frac_trros_k_500ns = calculate_hsasa_fraction(hsasa_trros_k_500ns, sasa_trros_k_500ns)


hsasa_frac_all_wt_500ns = calculate_hsasa_fraction(hsasa_all_wt_500ns, sasa_all_wt_500ns)

hsasa_frac_all_k_500ns = calculate_hsasa_fraction(hsasa_all_k_500ns, sasa_all_k_500ns)

# %%============================================================================
#                   Plot Rg-SASA KDE+cumulative
# =============================================================================

plot_kde_cumulative_seaborn(hsasa_frac_all_wt_500ns, hsasa_frac_all_k_500ns, figname='hsasa_sasa_kde_generated_encodermap',
                            type='SASA_nocumsum', name='')

# %%============================================================================
#                       Retrieve KDE data
# =============================================================================

df_wt = pd.DataFrame({'SASA': hsasa_frac_all_wt_500ns})
df_k = pd.DataFrame({'SASA': hsasa_frac_all_k_500ns})

df_wt_k = pd.concat([df_wt, df_k], ignore_index=True)
df_wt_k['Variant'] = np.repeat(['wt', 'K141E'], int(len(gyr_all_wt_500ns)))

sns.set_style('darkgrid')
fig, ax = plt.subplots(1, 1, figsize=(50, 20), gridspec_kw=dict(width_ratios=[4]))

ax = sns.kdeplot(data=df_wt_k, x="SASA", hue="Variant", clip=(0.0, 180.0), bw_adjust=0.8)

x = ax.lines[1].get_xdata()  # Get the x data of the distribution
y = ax.lines[1].get_ydata()  # Get the y data of the distribution
peaks, _ = np.array(find_peaks(y))

for p in range(len(peaks)):
    if y[peaks[p]] >= 0.015:
        print(x[peaks[p]])


# %%============================================================================
#                           Cumsum Rg
# =============================================================================

# cumsum_plot([gyr_trros_wt_500ns, gyr_trros_k_500ns], "Rg Cumulative frequency", "Rg (nm)",
#             True, "gyration_cumsum_wt_k141e")
#
# print(np.sum([gyr_trros_wt_500ns < 2.1]) / len(gyr_trros_wt_500ns))
# print(np.sum([gyr_trros_k_500ns < 2.1]) / len(gyr_trros_k_500ns))

# %%============================================================================
#                           Peaks Rg
# =============================================================================

# from scipy.signal import find_peaks
#
# sns.set_style('darkgrid')
#
# fig = plt.figure(figsize=(50, 40))
# ax = fig.add_subplot(111)
# plt.ylim(0, 0.05)
#
# (counts_hist_wt, hbins_wt, patches_wt) = ax.hist(gyr_all_wt_500ns, 100, density=True, stacked=True, alpha=0.6)
# for item in patches_wt:
#     item.set_height(item.get_height() / sum(counts_hist_wt))
#
# fig = plt.figure(figsize=(50, 40))
# ax = fig.add_subplot(111)
# plt.ylim(0, 0.05)
#
# (counts_hist_k, hbins_k, patches_k) = ax.hist(gyr_all_k_500ns, 100, density=True, stacked=True, color='orange',
#                                               alpha=0.6)
# for item in patches_k:
#     item.set_height(item.get_height() / sum(counts_hist_k))
#
# peaks_wt, _ = np.array(find_peaks(counts_hist_wt))
#
# peaks_k, _ = np.array(find_peaks(counts_hist_k))
#
# wt = [6, 13, 23, 34, 39, 41, 46, 49, 51, 56, 59, 62, 73, 75, 77, 81, 84, 86, 89, 92, 97]
# k = [3, 8, 12, 17, 23, 28, 30, 36, 38, 40, 45, 47, 49, 52, 59, 62, 64, 69, 75, 77, 80, 85, 89, 92, 95]
#
# num_frames_all = 15001
# frames_all_wt = np.linspace(1, num_frames_all, num_frames_all)
# frames_all_k = frames_all_wt
#
# np.set_printoptions(suppress=True)
# np.savetxt('gyration_peaks_delimiters_all_k.csv', [hbins_k[peaks_k - 1], hbins_k[peaks_k + 1]], delimiter=",", fmt='%f')
#
# for i in range(len(k)):
#     print(np.mean([hbins_k[peaks_k[i] - 1], hbins_k[peaks_k[i] + 1]]))
#     frames_peaks_all_k = [
#         frames_all_k[(gyr_all_k_500ns >= hbins_k[peaks_k[i] - 1]) & (gyr_all_k_500ns <= hbins_k[peaks_k[i] + 1])]]
#     print(frames_peaks_all_k)
#
#     pd.DataFrame(frames_peaks_all_k[0]).to_csv("gyration_peaks_all_k_{}.ndx".format(i), header=None, index=None)
#
# %% =============================================================================
#                                   Peaks Rg KDE
# =============================================================================

# peak1_rg_wt = df_wt.loc[(df_wt['Rg'] > 1.8895) & (df_wt['Rg'] < 1.8905)].index.to_list()
# peak2_rg_wt = df_wt.loc[(df_wt['Rg'] > 2.0895) & (df_wt['Rg'] < 2.0905)].index.to_list()
# peak3_rg_wt = df_wt.loc[(df_wt['Rg'] > 2.2995) & (df_wt['Rg'] < 2.3005)].index.to_list()
#
# peak1_rg_k = df_k.loc[(df_k['Rg'] > 1.8895) & (df_k['Rg'] < 1.8950)].index.to_list()
# peak2_rg_k = df_k.loc[(df_k['Rg'] > 2.0695) & (df_k['Rg'] < 2.0705)].index.to_list()
# peak3_rg_k = df_k.loc[(df_k['Rg'] > 2.6195) & (df_k['Rg'] < 2.6205)].index.to_list()
# peak4_rg_k = df_k.loc[(df_k['Rg'] > 2.6895) & (df_k['Rg'] < 2.6905)].index.to_list()
#
# print(peak4_rg_k)
