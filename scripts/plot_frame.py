import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Circle

fig, ax = plt.subplots()
ax.set(xlim=(-10, 110), ylim=(-10, 110))



#function to plot a circle
def create_sphere(x, y):

    coordinates = np.array([x, y])
    circle = Circle(xy = coordinates, radius = 1)
    circle.set_facecolor('red')
    ax.add_patch(circle)

#function to plot velocities
def velocity(x0, y0, u, v):

    #ax.quiver(x0, y0, u, v, color='b')
    x = [x0, u]
    y = [y0, v]
    ax.plot(x,y)

#plot forces
def forces(x0, y0, fx, fy):

    scale = 1E7
    ax.arrow(x0, y0, -fx/scale, -fy/scale, head_width=1, head_length=2)

#import data

data = pd.read_csv('results/frame_0.csv', sep=',',header=None, skiprows=1)

cols = len(data.columns)
rows = len(data)

print(rows)

a = 1

scale = 4.854368932E9

for i in range(0, rows):

    create_sphere(scale*data.loc[i, 0], scale*data.loc[i, 1])
    #velocity(data.loc[i, 0], data.loc[i, 1], data.loc[i, 3], data.loc[i, 4])
    #forces(data.loc[i, 0], data.loc[i, 1], data.loc[i, 6], data.loc[i, 7])



plt.show()
