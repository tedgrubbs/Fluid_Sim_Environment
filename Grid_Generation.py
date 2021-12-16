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
    STATIONARY_MOMENTUM_BASED = 5

# Use this to quickly redefine grid and config variables

D = 20

grid_size_x = 35*D+2
grid_size_y = 5*D+2
rho = np.zeros((grid_size_x,grid_size_y))
u = np.zeros((grid_size_x,grid_size_y))
v = np.zeros((grid_size_x,grid_size_y))
region = np.zeros((grid_size_x,grid_size_y),dtype='int')
boundary_v = np.zeros((grid_size_x,grid_size_y))
indices = np.zeros((grid_size_x,grid_size_y,2),dtype='int')
for i in range(grid_size_x):
    for j in range(grid_size_y):
        indices[i,j] = [i,j]

rho[:,:] = 1.
u[:,:] = 1.0

# creates border for sim environment
# top wall
region[:,-1] = Region.EXTERNAL
rho[:,-1] = 0.

# bottom wall
region[:,0] = Region.EXTERNAL
rho[:,0] = 0.

# right wall
region[-1, :] = Region.EXTERNAL
rho[-1, :] = 0.


# left wall
region[0, :] = Region.EXTERNAL
rho[0, :] = 0.

# left moving lid
# region[1, 1:-1] = 2
# region[2:-1,1] = 1
# region[2:-1,-2] = 1
# region[-2, 1:-1] = 1

# right moving lid
# region[-2, 1:-1] = 2
# region[1:-2,1] = 1
# region[1:-2,-2] = 1
# region[1, 1:-1] = 1

# bottom moving lid
# region[-2, 1:-1] = 1
# region[1:-2,-2] = 1
# region[1, 1:-1] = 1
# region[1:-1,1] = 2

# top moving lid with in and outlet
# region[-2, 1:-1] = Region.OUTLET
# region[1, 1:-1] = Region.INLET
# region[1:-2,1] = Region.STATIONARY
# region[1:-2,-2] = Region.MOVING_LID

# top moving lid with in and outlet. This configuration reproduce Borg's result when using only forward differences with predictor step
# region[-2, 1:-2] = Region.STATIONARY_MOMENTUM_BASED
# region[1, 1:-2] = Region.STATIONARY_MOMENTUM_BASED
# region[2:-2,1] = Region.MOVING_LID
# region[1:-1,-2] = Region.STATIONARY_MOMENTUM_BASED
# boundary_v[2:-2,1] = 1.0

# box in center
centerx = 5*D + D//2
centery = grid_size_y // 2

region[centerx-D//2 : centerx+D//2, centery-D//2 : centery+D//2] = Region.EXTERNAL
rho[centerx-D//2 : centerx+D//2, centery-D//2 : centery+D//2] = 0.
#
region[centerx-D//2 , centery-D//2-1 : centery+D//2+1] = Region.STATIONARY_MOMENTUM_BASED
region[centerx-D//2 : centerx+D//2+1, centery+D//2] = Region.STATIONARY_MOMENTUM_BASED
region[centerx+D//2 ,centery-D//2-1 : centery+D//2] = Region.STATIONARY_MOMENTUM_BASED
region[centerx-D//2 : centerx+D//2, centery-D//2-1] = Region.STATIONARY_MOMENTUM_BASED

region[-2, 1:-1] = Region.OUTLET

region[1, 2:-2] = Region.INLET
boundary_v[1, 2:-2] = 1.0

region[1:-2,1] = Region.MOVING_LID
boundary_v[1:-2,1] = 1.0

region[1:-2,-2] = Region.MOVING_LID
boundary_v[1:-2,-2] = 1.0

fig = plt.figure()
ax = fig.add_subplot(111)

ax.imshow(region.transpose())
ax.invert_yaxis()
plt.show()
# Initializing density everywhere

# u[:,-2] = 1.




output = pd.DataFrame(columns=['xi','yi','rho','u','v','region','boundary_v'])

output['xi'] = indices.reshape((-1,2))[:,0]
output['yi'] = indices.reshape((-1,2))[:,1]
output['rho'] = rho.reshape(-1)
output['u'] = u.reshape(-1)
output['v'] = v.reshape(-1)
output['region'] = region.reshape(-1)
output['boundary_v'] = boundary_v.reshape(-1)

output.to_csv('grid_variables.csv',index=False)


# Config variables
config = {}
config['grid_size_x'] = grid_size_x
config['grid_size_y'] = grid_size_y
config['real_size_x'] = (grid_size_x-2) // (grid_size_y-2)
config['real_size_y'] = 1.
config['frame_rate'] = 0
config['dt'] = 0.000175
config['dx'] = 1./(grid_size_x-3)*config['real_size_x']
config['dy'] = 1./(grid_size_y-3)*config['real_size_y']
config['viscosity'] = 1.8e-5
config['c'] = 347.0
config['force'] = 0.
config['run_graphics'] = 1
config['render_grid_size_x'] = 512*config['real_size_x']
config['render_grid_size_y'] = 256
config["tolerance"] = 0.00
config["max_run_time"] = 100000
config['thread_count'] = 4
config['Reynolds'] = 1000.
config['Mach'] = 0.1

with open('config.json','w') as fp:
    json.dump(config, fp, indent='\t')
