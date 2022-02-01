import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# %%============================================================================
#                           Projected DataFrame
# =============================================================================

projected = np.load("projected.npy")

projected_df = pd.DataFrame({'x': projected[0:, 0], 'y': projected[0:, 1]})
# projected_df['Variant'] = np.repeat(['wt'], int(len(projected[:,0])/2))

projected_df['Variant'] = np.repeat(['wt', 'K141E'], int(len(projected[:, 0]) / 2))

sns.set_style('darkgrid')

fig, ax = plt.subplots()

sns.kdeplot(data=projected_df, x='x', y='y', hue='Variant', fill=True, bw_adjust=0.2, legend=False)

ax.set_title("", fontsize=2)
ax.set_xlabel('x', fontsize=20)
ax.set_ylabel('y', fontsize=20)

for label in (ax.get_xticklabels() + ax.get_yticklabels()):
    label.set_fontsize(16)

# GENERATED ENCODERMAP STRUCTURES
plt.scatter(-3.748, 0.850, marker="x", c="r")
plt.scatter(-4.494, -1.311, marker="x", c="r")
plt.scatter(-3.419, -0.976, marker="x", c="r")
plt.scatter(-4.011, -0.049, marker="x", c="r")
plt.scatter(-2.008, -0.364, marker="x", c="r")
plt.scatter(-2.484, -2.335, marker="x", c="r")
plt.scatter(-0.985, -1.828, marker="x", c="r")
plt.scatter(-0.247, -0.584, marker="x", c="r")
plt.scatter(0.419, -0.020, marker="x", c="r")
plt.scatter(0.126, -1.101, marker="x", c="r")
plt.scatter(-1.577, 1.109, marker="x", c="r")
plt.scatter(-1.277, 1.558, marker="x", c="r")
plt.scatter(-0.649, 1.606, marker="x", c="r")
plt.scatter(-0.371, 1.166, marker="x", c="r")
plt.scatter(0.638, 2.706, marker="x", c="r")
plt.scatter(1.413, 3.500, marker="x", c="r")
plt.scatter(2.283, 3.127, marker="x", c="r")
plt.scatter(1.479, 1.749, marker="x", c="r")
plt.scatter(2.619, 0.296, marker="x", c="r")
plt.scatter(4.337, 3.414, marker="x", c="r")
plt.scatter(4.681, 2.438, marker="x", c="r")
plt.scatter(5.134, 2.400, marker="x", c="r")

plt.show()

fig.savefig("KDE_encodermap_crosses_paper.png", dpi=600)
