import numpy as np
import pandas as pd
import json
from matplotlib import pyplot as plt
plt.style.use('dark_background')

# Use this to quickly redefine grid and config variables

grid_size_x = 256
grid_size_y = 256
rho = np.zeros((grid_size_x,grid_size_y))
u = np.zeros((grid_size_x,grid_size_y))
v = np.zeros((grid_size_x,grid_size_y))
boundary = np.zeros((grid_size_x,grid_size_y),dtype='int')
indices = np.zeros((grid_size_x,grid_size_y,2),dtype='int')
for i in range(grid_size_x):
    for j in range(grid_size_y):
        indices[i,j] = [i,j]

# Top and bottom walls
boundary[:,0] = 1
boundary[:,-1] = 2

# left and right walls
boundary[0,:-1] = 1
boundary[-1,:-1] = 1

# Initializing density everywhere
rho[:,:] = 1.
u[:,-1] = 1.

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
config['dx'] = 1./grid_size_x*config['real_size_x']
config['dy'] = 1./grid_size_y*config['real_size_x']
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
