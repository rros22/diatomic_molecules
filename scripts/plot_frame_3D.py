import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.patches import Circle


fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

ax.set_xlim3d(0, 100)
ax.set_ylim3d(0, 100)
ax.set_zlim3d(0, 100)


def velocity(x0, y0, z0, u, v, w):

    ax.quiver(x0, y0, z0, u, v, w)


#function to plot a circle
def create_sphere(x0, y0, z0, r):

    u = np.linspace(0, 2 * np.pi, 100)
    v = np.linspace(0, np.pi, 100)
    x = r * np.outer(np.cos(u), np.sin(v)) + x0
    y = r * np.outer(np.sin(u), np.sin(v)) + y0
    z = r * np.outer(np.ones(np.size(u)), np.cos(v)) + z0

    ax.plot_surface(x, y, z, color='b')




#import data

data = pd.read_csv('results/frame.csv', sep=',',header=None, skiprows=1)

cols = len(data.columns)
rows = len(data)

for i in range(0, rows):
    create_sphere(data.loc[i, 0], data.loc[i, 1], data.loc[i, 2], 2)
    velocity(data.loc[i, 0], data.loc[i, 1], data.loc[i, 2], data.loc[i, 3], data.loc[i, 4], data.loc[i, 5])


plt.show()
