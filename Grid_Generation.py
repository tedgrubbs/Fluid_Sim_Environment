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
    STATIC_U = 6
    EX_OUTLET = 7

# Use this to quickly redefine grid and config variables

D = 70

# grid_size_x = 35*D+2
# grid_size_y = 3*D+2
grid_size_x = D+2
grid_size_y = D+2

rho = np.zeros((grid_size_x,grid_size_y))
u = np.zeros((grid_size_x,grid_size_y))
v = np.zeros((grid_size_x,grid_size_y))
region = np.zeros((grid_size_x,grid_size_y),dtype='int')
boundary_v = np.zeros((grid_size_x,grid_size_y))
indices = np.zeros((grid_size_x,grid_size_y,2),dtype='int')
for i in range(grid_size_x):
    for j in range(grid_size_y):
        indices[i,j] = [i,j]

rho[:,:] = 1.22
# u[:,:] = 1.0

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

# This configuration reproduce Borg's result when using only forward differences with predictor step and using basic STATIONARY region type
# Note that at the corners where the moving lid intersects the stationary walls, these should be marked as Moving lid points.
# Otherwise the density at the corners will grow indefinitely- even though the rest of the simulation is stable. This is how Borg's simulations works. 
# This is caused by this term in the density equation: (-ru[i-2][j] + 4.*ru[i-1][j] - 3.*ru[i][j])
region[-2, 1:-1] = Region.STATIONARY_MOMENTUM_BASED
region[1, 1:-1] = Region.STATIONARY_MOMENTUM_BASED
region[1:-1,-2] = Region.MOVING_LID
region[1:-1,1] = Region.STATIONARY_MOMENTUM_BASED
boundary_v[1:-1,-2] = 34.7

# u[1:-1,-2] = 1.0

# flow over flat plate. Be sure to turn down timestep for this at high mach number
# region[-2, 1:-1] = Region.EX_OUTLET
# region[1, 1:-2] = Region.STATIC_U
# region[1:-2,-2] = Region.STATIC_U
# region[1:-2,1] = Region.STATIONARY_MOMENTUM_BASED
# boundary_v[1, 2:-2] = 1.0
# boundary_v[1:-2,-2] = 1.0
# u[1:-1,1:] = 1.0

# box in center
# centerx = 15*D + D//2
# centery = grid_size_y // 2

# region[centerx-D//2 : centerx+D//2, centery-D//2 : centery+D//2] = Region.EXTERNAL

# region[centerx-D//2 , centery-D//2-1 : centery+D//2+1] = Region.STATIONARY_MOMENTUM_BASED
# region[centerx-D//2 : centerx+D//2+1, centery+D//2] = Region.STATIONARY_MOMENTUM_BASED
# region[centerx+D//2 ,centery-D//2-1 : centery+D//2] = Region.STATIONARY_MOMENTUM_BASED
# region[centerx-D//2 : centerx+D//2, centery-D//2-1] = Region.STATIONARY_MOMENTUM_BASED

# region[-2, 1:-2] = Region.OUTLET
# boundary_v[-2, 2:-2] = 1.0

# region[1, 2:-2] = Region.INLET
# boundary_v[1, 2:-2] = 1.0

# region[1:-1,1] = Region.MOVING_LID
# boundary_v[1:-1,1] = 1.0

# region[1:-1,-2] = Region.MOVING_LID
# boundary_v[1:-1,-2] = 1.0

# side block
# region[centerx-D//2+1 : centerx+D//2, centery+D//2 : -1] = Region.EXTERNAL
# region[centerx-D//2 , centery+D//2 : -2] = Region.STATIONARY_MOMENTUM_BASED
# region[centerx-D//2 : centerx+D//2+1, centery+D//2] = Region.STATIONARY_MOMENTUM_BASED
# region[centerx+D//2 , centery+D//2 : -2] = Region.STATIONARY_MOMENTUM_BASED
# region[centerx-D//2 : centerx+D//2, centery-D//2-1] = Region.STATIONARY_MOMENTUM_BASED


fig = plt.figure()
ax = fig.add_subplot(111)

ax.imshow(region.transpose())
ax.invert_yaxis()
plt.show()




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
config['real_size_y'] = 0.01
config['real_size_x'] = 0.01



config['frame_rate'] = 0
config['dt'] = 0.000175
config['dx'] = 1./(grid_size_x-3)*config['real_size_x']
config['dy'] = 1./(grid_size_y-3)*config['real_size_y']
config['viscosity'] = 1.8e-5
config['c'] = 347.0
config['force'] = 0.
config['run_graphics'] = 1

base_render = 512
if grid_size_x > base_render:
    base_render = grid_size_x
    base_render = np.clip(base_render, a_min=None, a_max=1024)
y_multiplier = config['real_size_y'] / config['real_size_x']
config['render_grid_size_x'] = int(base_render)
config['render_grid_size_y'] = int(base_render*y_multiplier)

config["tolerance"] = 0.00
config["max_run_time"] = 100000
config['thread_count'] = 4
config['Reynolds'] = 400.
config['Mach'] = 0.1

print('Reynolds number:', rho[1,-2]*boundary_v[1,-2]*config['real_size_x']/config['viscosity'])

with open('config.json','w') as fp:
    json.dump(config, fp, indent='\t')
