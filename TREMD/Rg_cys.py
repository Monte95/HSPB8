from remd_distributions import *

t, gyr_all_wt = read_gyration_data("gyrate_all_wt.xvg")
t, gyr_all_k = read_gyration_data("gyrate_all_k.xvg")

t, gyr_all_wt_10_99 = read_gyration_data("gyrate_wt_10_99.xvg")
t, gyr_all_wt_10_195 = read_gyration_data("gyrate_wt_10_195.xvg")
t, gyr_all_wt_99_195 = read_gyration_data("gyrate_wt_99_195.xvg")

t, gyr_all_k_10_99 = read_gyration_data("gyrate_k_10_99.xvg")
t, gyr_all_k_10_195 = read_gyration_data("gyrate_k_10_195.xvg")
t, gyr_all_k_99_195 = read_gyration_data("gyrate_k_99_195.xvg")

plot_kde_cumulative_seaborn(gyr_all_wt_99_195, gyr_all_k_99_195, figname='gyration_kde_wt_k_99_195',
                            type='Rg_nocumsum', name='')

plt.show()

