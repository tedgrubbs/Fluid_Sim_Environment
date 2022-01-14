import numpy as np
import pandas as pd
import json
from matplotlib import pyplot as plt
from enum import IntEnum
plt.style.use('dark_background')

class Region(IntEnum):
    EXTERNAL = -1
    FREE_FLOW = 0
    LEFT_WALL = 1
    RIGHT_WALL = 2
    TOP_WALL = 3
    BOTTOM_WALL = 4
    TOP_MOVING_LID = 5
    STATIC = 6
    RIGHT_OUTFLOW = 7
    LEFT_INLET = 8
    RIGHT_PRESSURE_OUTLET = 9
    CORNER_POINT = 10
    PERIODIC_Y_TOP = 11
    PERIODIC_Y_BOTTOM = 12
    

# Use this to quickly redefine grid and config variables

D = 30
SPEED = 340.28 # speed of sound at STP

# grid_size_x = 35*D+2
# grid_size_y = 3*D+2

grid_size_x = D*30+2
grid_size_y = D+2

# grid_size_x = D+2
# grid_size_y = D+2

rho = np.zeros((grid_size_x,grid_size_y))
u = np.zeros((grid_size_x,grid_size_y))
v = np.zeros((grid_size_x,grid_size_y))
temperature = np.zeros((grid_size_x,grid_size_y))
region = np.zeros((grid_size_x,grid_size_y),dtype='int')
indices = np.zeros((grid_size_x,grid_size_y,2),dtype='int')
for i in range(grid_size_x):
    for j in range(grid_size_y):
        indices[i,j] = [i,j]

rho[:,:] = 1.22
temperature[:,:] = 288
# temperature[int(D/2)-5:int(D/2)+5, int(D/2)-5:int(D/2)+5] = 310.
# u[:,:] = 1.0

# creates border for sim environment
# top wall
region[:,-1] = Region.EXTERNAL

# bottom wall
region[:,0] = Region.EXTERNAL

# right wall
region[-1, :] = Region.EXTERNAL

# left wall
region[0, :] = Region.EXTERNAL

# This configuration reproduce Borg's result when using only forward differences with predictor step and using basic STATIONARY region type
# Note that at the corners where the moving lid intersects the stationary walls, these should be marked as Moving lid points.
# Otherwise the density at the corners will grow indefinitely- even though the rest of the simulation is stable. This is how Borg's simulations works. 
# This is caused by this term in the density equation: (-ru[i-2][j] + 4.*ru[i-1][j] - 3.*ru[i][j])
# region[-2, 1:-1] = Region.RIGHT_WALL
# region[1, 1:-2] = Region.LEFT_WALL
# region[1:-1,-2] = Region.TOP_MOVING_LID
# region[1:-2,1] = Region.BOTTOM_WALL
# u[1:-1,-2] = 0.1*SPEED

# region[-2, 1:-1] = Region.RIGHT_WALL
# region[1, 1:-1] = Region.LEFT_WALL
# region[2:-2,-2] = Region.PERIODIC_Y_TOP
# region[2:-2,1] = Region.PERIODIC_Y_BOTTOM
# rho[int(D/2)-5:int(D/2)+5, int(10)-5:int(10)+5] = 2.0

centerx = int(0.222222222 * grid_size_x)
centery = int(2./3.* grid_size_y)
length = int(0.055555556 * grid_size_x)

# Left Velocity inlet, right outflow
region[-2, 1:-1] = Region.RIGHT_PRESSURE_OUTLET
region[1, 2:-2] = Region.LEFT_INLET
region[1:-2,-2] = Region.PERIODIC_Y_TOP
region[1:-2,1] = Region.PERIODIC_Y_BOTTOM

# Creating a box in the flow path
region[centerx, 2 : centery] = Region.RIGHT_WALL

region[centerx+1:centerx+length, 1 : centery] = Region.EXTERNAL

region[centerx+length, 2 : centery] = Region.LEFT_WALL

region[centerx+1:centerx+length , centery] = Region.BOTTOM_WALL

region[centerx:centerx+length+1, -2] = Region.TOP_WALL

region[centerx , centery] = Region.CORNER_POINT
region[centerx+length , centery] = Region.CORNER_POINT

u[1, 2:-2] = .1

# flow over flat plate. Be sure to turn down timestep for this at high mach number
# region[1:-2,1] = Region.BOTTOM_WALL
# region[-2, 1:-1] = Region.RIGHT_OUTFLOW
# region[1, 2:-2] = Region.STATIC
# region[1:-2,-2] = Region.STATIC
# u[1, 2:-2] = SPEED * 4.
# u[1:-2,-2] = SPEED * 4.
# u[1:-1,1:-1] = SPEED * 4.
# u[1:-2,1] = 0.

# box in center
# centerx = 15*D + D//2
# centery = grid_size_y // 2

# region[centerx-D//2 : centerx+D//2, centery-D//2 : centery+D//2] = Region.EXTERNAL

# region[centerx-D//2 , centery-D//2-1 : centery+D//2+1] = Region.STATIONARY_MOMENTUM_BASED
# region[centerx-D//2 : centerx+D//2+1, centery+D//2] = Region.STATIONARY_MOMENTUM_BASED
# region[centerx+D//2 ,centery-D//2-1 : centery+D//2] = Region.STATIONARY_MOMENTUM_BASED
# region[centerx-D//2 : centerx+D//2, centery-D//2-1] = Region.STATIONARY_MOMENTUM_BASED

# region[-2, 1:-2] = Region.OUTLET

# region[1, 2:-2] = Region.INLET

# region[1:-1,1] = Region.MOVING_LID

# region[1:-1,-2] = Region.MOVING_LID

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




output = pd.DataFrame(columns=['xi','yi','rho','u','v','temperature','region'])

output['xi'] = indices.reshape((-1,2))[:,0]
output['yi'] = indices.reshape((-1,2))[:,1]
output['rho'] = rho.reshape(-1)
output['u'] = u.reshape(-1)
output['v'] = v.reshape(-1)
output['temperature'] = temperature.reshape(-1)
output['region'] = region.reshape(-1)

output.to_csv('grid_variables.csv',index=False)


# Config variables
config = {}
config['grid_size_x'] = grid_size_x
config['grid_size_y'] = grid_size_y
config['real_size_y'] = 0.03#/1.218487395
config['real_size_x'] = 0.9



config['frame_rate'] = 0
config['dt'] = 0.000175
config['dx'] = 1./(grid_size_x-3)*config['real_size_x']
config['dy'] = 1./(grid_size_y-3)*config['real_size_y']
config['viscosity'] = 1.81e-5
config['c'] = SPEED
config['gamma'] = 1.4
config['Pr'] = 0.071
config['run_graphics'] = 1

base_render = 512
if grid_size_x > base_render:
    base_render = grid_size_x
    base_render = np.clip(base_render, a_min=None, a_max=1024)
y_multiplier = config['real_size_y'] / config['real_size_x']
config['render_grid_size_x'] = int(base_render)
config['render_grid_size_y'] = int(base_render*y_multiplier)*5

config["tolerance"] = 0.00
config["max_run_time"] = 2000000
config['thread_count'] = 4

print('Reynolds number:', rho[1,-2]*u[1,2]*config['real_size_y']/config['viscosity'])
print('C Reynolds number:', rho[1,-2]*SPEED*config['real_size_x']/config['viscosity'])
# print('Grid x Reynolds number:', rho[1,2]*u[1, 2]*config['dx']/config['viscosity'])
# print('Grid y Reynolds number:', rho[1,2]*u[1, 2]*config['dy']/config['viscosity'])

with open('config.json','w') as fp:
    json.dump(config, fp, indent='\t')
