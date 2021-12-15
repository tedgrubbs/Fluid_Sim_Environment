import numpy as np
import pandas as pd
import json
from matplotlib import pyplot as plt
from enum import IntEnum
plt.style.use('dark_background')

class Region(IntEnum):
    EXTERNAL = -1
    FREE_FLOW = 0
    STATIONARY = 1
    MOVING_LID = 2
    INLET = 3
    OUTLET = 4

# Use this to quickly redefine grid and config variables

grid_size_x = 65+2
grid_size_y = 65+2
rho = np.zeros((grid_size_x,grid_size_y))
u = np.zeros((grid_size_x,grid_size_y))
v = np.zeros((grid_size_x,grid_size_y))
boundary = np.zeros((grid_size_x,grid_size_y),dtype='int')
indices = np.zeros((grid_size_x,grid_size_y,2),dtype='int')
for i in range(grid_size_x):
    for j in range(grid_size_y):
        indices[i,j] = [i,j]

rho[:,:] = 1.

# top wall
boundary[:,-1] = Region.EXTERNAL
rho[:,-1] = 0.

# bottom wall
boundary[:,0] = Region.EXTERNAL
rho[:,0] = 0.

# right wall
boundary[-1, :] = Region.EXTERNAL
rho[-1, :] = 0.


# left wall
boundary[0, :] = Region.EXTERNAL
rho[0, :] = 0.

# left moving lid
# boundary[1, 1:-1] = 2
# boundary[2:-1,1] = 1
# boundary[2:-1,-2] = 1
# boundary[-2, 1:-1] = 1

# right moving lid
# boundary[-2, 1:-1] = 2
# boundary[1:-2,1] = 1
# boundary[1:-2,-2] = 1
# boundary[1, 1:-1] = 1

# bottom moving lid
# boundary[-2, 1:-1] = 1
# boundary[1:-2,-2] = 1
# boundary[1, 1:-1] = 1
# boundary[1:-1,1] = 2

# top moving lid
boundary[-2, 1:-1] = Region.STATIONARY
boundary[1, 1:-1] = Region.STATIONARY
boundary[1:-1,1] = Region.STATIONARY
boundary[2:-2,-2] = Region.MOVING_LID

# left and right moving lid
# boundary[1, 1:-1] = 2
# boundary[-2, 1:-1] = 2
# boundary[2:-2,1] = 1
# boundary[2:-2,-2] = 1

# all lids moving
# boundary[1, 1:-1] = 2
# boundary[-2, 1:-1] = 2
# boundary[2:-2,1] = 2
# boundary[2:-2,-2] = 2

# 3 lids moving
# boundary[1, 1:-1] = 2
# boundary[-2, 1:-1] = 2
# boundary[2:-2,1] = 1
# boundary[2:-2,-2] = 2

# top lid moving with barrier in center
# boundary[1:-1,-2] = 1
# boundary[1, 1:-2] = 1
# boundary[-2, 1:-2] = 1
# boundary[1:-2, 1] = 1

# boundary[33, :] = -1
# boundary[34, 1:-1] = 2
# boundary[32, 1:-1] = 2

fig = plt.figure()
ax = fig.add_subplot(111)

ax.imshow(boundary.transpose())
ax.invert_yaxis()
plt.show()
# Initializing density everywhere

# u[:,-2] = 1.




output = pd.DataFrame(columns=['xi','yi','rho','u','v','boundary'])

output['xi'] = indices.reshape((-1,2))[:,0]
output['yi'] = indices.reshape((-1,2))[:,1]
output['rho'] = rho.reshape(-1)
output['u'] = u.reshape(-1)
output['v'] = v.reshape(-1)
output['boundary'] = boundary.reshape(-1)

output.to_csv('grid_variables.csv',index=False)


# Config variables
config = {}
config['grid_size_x'] = grid_size_x
config['grid_size_y'] = grid_size_y
config['real_size_x'] = 1.
config['real_size_y'] = 1.
config['frame_rate'] = 0
config['dt'] = 0.000175
config['dx'] = 1./(grid_size_x-3)*config['real_size_x']
config['dy'] = 1./(grid_size_y-3)*config['real_size_y']
config['viscosity'] = 1.8e-5
config['c'] = 347.0
config['force'] = 0.
config['run_graphics'] = 1
config['render_grid_size_x'] = 512
config['render_grid_size_y'] = 512
config["tolerance"] = 0.00
config["max_run_time"] = 100000

with open('config.json','w') as fp:
    json.dump(config, fp, indent='\t')
