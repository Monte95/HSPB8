import numpy as np
import matplotlib.pyplot as plt

projected = np.load('projected.npy')

fig1, axe1 = plt.subplots()

hist, xedges, yedges = np.histogram2d(projected[:, 0], projected[:, 1], bins=200)

caxe = axe1.imshow(-np.log(hist.T), origin='lower', extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]], aspect="auto")
cbar = fig1.colorbar(caxe)
cbar.set_label("-ln(p)", labelpad=0)
axe1.set_title("Path Generator")

plt.show()
np.set_printoptions(threshold=np.inf)

print(hist)

