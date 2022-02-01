import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import pandas as pd
from scipy import stats


def read_time_contacts_data(document_name):
    """Reads contacts in every frame between protein and paroxetine from contact_A.dat obtained from newcontact.tcl
    and creates a dataframe"""

    time, residue, distance = [], [], []
    with open(document_name) as f:
        for line in f:
            cols = line.split()

            if len(cols) == 8:
                time.append(int(int(cols[0]) / 10))
                residue.append(int(cols[1]))
                distance.append(float(cols[7]))

    time_contacts_data_df = pd.DataFrame({'Time': np.array(time),
                                          'Residue': np.array(residue),
                                          'Distance': np.array(distance)})

    return time_contacts_data_df


def read_frequency_contacts_data(document_name):
    """Reads contact frequency data from contact_A.dat obtained from newcontact.tcl and creates a dataframe ordered by
    residue number"""

    residue_num, residue_name, frequency = [], [], []
    with open(document_name) as f:
        for line in f:
            cols = line.split()

            if len(cols) == 7:
                residue_num.append(int(cols[0]) + 1)
                residue_name.append(cols[1])
                frequency.append(float(cols[3]))

    frequency_contacts_data_df = pd.DataFrame({'Residue Number': np.array(residue_num),
                                               'Residue Name': np.array(residue_name),
                                               'Frequency': np.array(frequency)})

    frequency_contacts_data_df = frequency_contacts_data_df.sort_values(by=['Residue Number'])

    return frequency_contacts_data_df


# CREATE CONTACT TIME DATAFRAMES
wt_paroxetine_contacts_time_df = read_time_contacts_data("contact_all_wt_paroxetine.dat")
k_paroxetine_contacts_time_df = read_time_contacts_data("contact_all_k_paroxetine.dat")

wt_k_paroxetine_contacts_time_df = pd.concat([wt_paroxetine_contacts_time_df,
                                              k_paroxetine_contacts_time_df], ignore_index=True)

wt_k_paroxetine_contacts_time_df['Variant'] = np.repeat(['wt', 'K141E'],
                                                        [len(wt_paroxetine_contacts_time_df['Distance']),
                                                         len(k_paroxetine_contacts_time_df['Distance'])])

# DOCKED TIME PERCENTAGE
wt_paroxetine_docked_time_percentage = len(wt_paroxetine_contacts_time_df['Time'].unique()) / \
                                          (wt_paroxetine_contacts_time_df['Time'].unique().max() -
                                           wt_paroxetine_contacts_time_df['Time'].unique().min() + 1) * 100

k_paroxetine_docked_time_percentage = len(k_paroxetine_contacts_time_df['Time'].unique()) / \
                                         (k_paroxetine_contacts_time_df['Time'].unique().max() -
                                          k_paroxetine_contacts_time_df['Time'].unique().min() + 1) * 100

docked_percentage_df = pd.DataFrame({'Docked': ['Docked', 'Not docked', 'Docked', 'Not docked'],
                                     'Percentage': [wt_paroxetine_docked_time_percentage,
                                                    100-wt_paroxetine_docked_time_percentage,
                                                    k_paroxetine_docked_time_percentage,
                                                    100-k_paroxetine_docked_time_percentage],
                                     'Variant': ['wt', 'wt', 'K141E', 'K141E']})

# CREATE CONTACT FREQUENCY DATAFRAMES
wt_paroxetine_contacts_frequency_df = read_frequency_contacts_data('contact_A_wt_paroxetine.dat')
k_paroxetine_contacts_frequency_df = read_frequency_contacts_data('contact_A_k_paroxetine.dat')

wt_k_paroxetine_contacts_frequency_df = pd.concat([wt_paroxetine_contacts_frequency_df,
                                                   k_paroxetine_contacts_frequency_df], ignore_index=True)

wt_k_paroxetine_contacts_frequency_df['Variant'] = np.repeat(['wt', 'K141E'],
                                                             [len(wt_paroxetine_contacts_frequency_df['Frequency']),
                                                              len(k_paroxetine_contacts_frequency_df['Frequency'])])

# PLOT DOCKED PERCENTAGES
sns.set_style('darkgrid')

g = sns.catplot(data=docked_percentage_df, x='Docked', y='Percentage', hue='Variant',
                kind='bar', legend_out=False, legend=False)

g.set(ylim=(0, 100))

ax = g.facet_axis(0, 0)

for c in ax.containers:
    labels = [f'{(v.get_height()):.1f} %' for v in c]
    ax.bar_label(c, labels=labels, label_type='edge', fontsize=16)

for label in (ax.get_xticklabels() + ax.get_yticklabels()):
    label.set_fontsize(16)

ax.set_xlabel('', size=16)
ax.set_ylabel('Percentage (%)', size=16)

ax.set_title("Paroxetine docking percentage", size=16)

leg = ax.legend(prop={"size": 18}, loc='upper right')

# plt.show()
plt.close()


# g.figure.savefig("paroxetine_docking_percentages_wt_k.png", dpi=320)


# PLOT TIMELINE CONTACT DATA SINGLE TRAJECTORIES
wt_7_paroxetine_contacts_time_df = wt_paroxetine_contacts_time_df.loc[wt_paroxetine_contacts_time_df['Time'] < 500]

g = sns.relplot(data=wt_7_paroxetine_contacts_time_df, x='Time', y='Residue', size=0.05, legend=False)

ticks = np.arange(1, 197, 1)
plt.yticks(ticks)

labels_left = np.arange(1, 197, 1).tolist()
labels_right = np.arange(1, 197, 1).tolist()

# Leave only even or odd labels
for i in range(0, len(labels_left)):
    if labels_left[i] % 2 != 0:
        labels_left[i] = 0

for i in range(0, len(labels_right)):
    if labels_right[i] % 2 == 0:
        labels_right[i] = 0

# Leave only labels that correspond to contacts
for i in range(0, len(labels_left)):
    if labels_left[i] not in wt_7_paroxetine_contacts_time_df['Residue'].to_list():
        labels_left[i] = ''

for i in range(0, len(labels_right)):
    if labels_right[i] not in wt_7_paroxetine_contacts_time_df['Residue'].to_list():
        labels_right[i] = ''

# Set new labels for left y-axis for even labels
g.set_yticklabels(labels_left, size=8)

# Set new labels for right y-axis for odd labels
ax2 = plt.twinx()
ax2.set_yticks(ticks)
ax2.set_yticklabels(labels_right, size=8)

# plt.show()
plt.close()


# PLOT TIMELINE CONTACT DATA
g = sns.relplot(data=wt_k_paroxetine_contacts_time_df, x='Time', y='Residue', size=0.05, row='Variant', hue='Variant',
                palette=["blue", "orange"], legend=False)

g.map(plt.axvline, x=500, color=".7", dashes=(2, 1), zorder=0)
g.map(plt.axvline, x=1000, color=".7", dashes=(2, 1), zorder=0)
g.map(plt.axvline, x=1500, color=".7", dashes=(2, 1), zorder=0)
g.map(plt.axvline, x=2000, color=".7", dashes=(2, 1), zorder=0)
g.map(plt.axvline, x=2500, color=".7", dashes=(2, 1), zorder=0)
g.map(plt.axvline, x=3000, color=".7", dashes=(2, 1), zorder=0)
g.map(plt.axvline, x=3500, color=".7", dashes=(2, 1), zorder=0)
g.map(plt.axvline, x=4000, color=".7", dashes=(2, 1), zorder=0)
g.map(plt.axvline, x=4500, color=".7", dashes=(2, 1), zorder=0)

# Titles and axis
axes = g.axes.flat
axes[0].set_title("wt", size=18)
axes[1].set_title("K141E", size=18)

axes[0].set_xlabel('', size=18)
axes[1].set_xlabel('Time (ns)', size=18)

plt.xticks(np.arange(0, 5000, 500))

axes[0].set_ylabel('Residue', size=18)
axes[1].set_ylabel('Residue', size=18)

# plt.show()
plt.close()


# g.figure.savefig("paroxetine_contact_time_wt_k.png", dpi=320)


# PLOT FREQUENCY DATA
colors = np.repeat(['Blue', 'Yellow', 'Red'], [95, 74, 27])

g = sns.catplot(data=wt_k_paroxetine_contacts_frequency_df, x='Residue Number', y='Frequency', row='Variant',
                kind='bar', legend_out=False, palette=colors, legend=False)

# Titles and axis
axes = g.axes.flat
axes[0].set_title("wt", size=18)
axes[1].set_title("K141E", size=18)

axes[0].set_xlabel('', size=18)
axes[1].set_xlabel('Residue', size=18)

axes[0].set_ylabel('Contact Frequency (%)', size=18)
axes[1].set_ylabel('Contact Frequency (%)', size=18)

# Reduce xticks frequency
for ax in g.axes.flat:

    labels = ax.get_xticklabels()  # get x labels

    for i, l in enumerate(labels):
        if i % 5 != 0:
            labels[i] = ''  # skip labels multiple of 5

    ax.set_xticklabels(labels, size=14)  # set new labels

# plt.show()

# g.figure.savefig("paroxetine_contact_frequency_wt_k.png", dpi=320)

plt.close()
