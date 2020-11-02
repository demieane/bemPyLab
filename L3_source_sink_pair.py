#--------------------------------------------------------------------------------
# TUTORIAL_3.1
# Visualize the velocity field of a source-sink pair (2d)
#--------------------------------------------------------------------------------

import numpy as np
from matplotlib import pyplot

N = 50                                # number of points in each direction
x_start, x_end = -2.0, 2.0            # boundaries in the x-direction
y_start, y_end = -1.0, 1.0            # boundaries in the y-direction
x = np.linspace(x_start, x_end, N)    # creates a 1D-array with the x-coordinates
y = np.linspace(y_start, y_end, N)    # creates a 1D-array with the y-coordinates
X, Y = np.meshgrid(x, y)              # generates a mesh grid

# plot the grid of points
#--------------------------------------------------------------------------------
width = 10.0
height = (y_end - y_start) / (x_end - x_start) * width
pyplot.figure(figsize=(width, height))
pyplot.xlabel('x', fontsize=16)
pyplot.ylabel('y', fontsize=16)
pyplot.xlim(x_start, x_end)
pyplot.ylim(y_start, y_end)
pyplot.title('meshgrid')
pyplot.scatter(X, Y, s=5, color='#CD2305', marker='o')
pyplot.show() #it tells the terminal to make the plot appear
#--------------------------------------------------------------------------------

strength_source = 5.0                      # source strength
x_source, y_source = -1.0, 0.0             # location of the source

# compute the velocity field on the mesh grid
u_source = (strength_source / (2 * np.pi) *
            (X - x_source) / ((X - x_source)**2 + (Y - y_source)**2))
v_source = (strength_source / (2 * np.pi) *
            (Y - y_source) / ((X - x_source)**2 + (Y - y_source)**2))

# plot the streamlines
#--------------------------------------------------------------------------------
width = 10.0
height = (y_end - y_start) / (x_end - x_start) * width
pyplot.figure(figsize=(width, height))
pyplot.xlabel('x', fontsize=16)
pyplot.ylabel('y', fontsize=16)
pyplot.xlim(x_start, x_end)
pyplot.ylim(y_start, y_end)
pyplot.title('streamlines - source')
pyplot.streamplot(X, Y, u_source, v_source, density=2, linewidth=1, arrowsize=2, arrowstyle='->')
pyplot.scatter(x_source, y_source, color='#CD2305', s=80, marker='o')
pyplot.show()
#--------------------------------------------------------------------------------

strength_sink = -5.0                     # strength of the sink
x_sink, y_sink = 1.0, 0.0                # location of the sink

# compute the velocity on the mesh grid
u_sink = (strength_sink / (2 * np.pi) *
          (X - x_sink) / ((X - x_sink)**2 + (Y - y_sink)**2))
v_sink = (strength_sink / (2 * np.pi) *
          (Y - y_sink) / ((X - x_sink)**2 + (Y - y_sink)**2))

# plot the streamlines
#--------------------------------------------------------------------------------
width = 10.0
height = (y_end - y_start) / (x_end - x_start) * width
pyplot.figure(figsize=(width, height))
pyplot.xlabel('x', fontsize=16)
pyplot.ylabel('y', fontsize=16)
pyplot.xlim(x_start, x_end)
pyplot.ylim(y_start, y_end)
pyplot.title('streamlines - sink')
pyplot.streamplot(X, Y, u_sink, v_sink,
                density=2, linewidth=1, arrowsize=2, arrowstyle='->')
pyplot.scatter(x_sink, y_sink,
            color='#CD2305', s=80, marker='o')
pyplot.show()
#--------------------------------------------------------------------------------

# compute the velocity of the pair source/sink by superposition
u_pair = u_source + u_sink
v_pair = v_source + v_sink

# plot the streamlines of the pair source/sink
#--------------------------------------------------------------------------------
width = 10.0
height = (y_end - y_start) / (x_end - x_start) * width
pyplot.figure(figsize=(width, height))
pyplot.xlabel('x', fontsize=16)
pyplot.ylabel('y', fontsize=16)
pyplot.xlim(x_start, x_end)
pyplot.ylim(y_start, y_end)
pyplot.title('streamlines - source/sink pair')
pyplot.streamplot(X, Y, u_pair, v_pair,
                density=2.0, linewidth=1, arrowsize=2, arrowstyle='->')
pyplot.scatter([x_source, x_sink], [y_source, y_sink], 
            color='#CD2305', s=80, marker='o')

pyplot.show()
#--------------------------------------------------------------------------------