import numpy as np

method = 'line'  # Factor or line
x_start = np.array([1.0, 1.0, 1.0, 1.0, 1.0, 1.0])  # Must be a numpy array


potential_name = 'Cu01.eam.alloy'
potential_type = 'eam/alloy'  # Will not work for other types yet
finite_size_offset = -252  # Difference between infinite system size and finite system size (infinite - finite)

# Must be the same number of prefixes and targets.
prefixes = ['32Cu-0.1', '32Cu-30']
targets = [1299, 2149]

# NM parameters
step = 0.05
no_improve_thr = 20
no_improv_break = 3
max_iter = 15
alpha = 1.
gamma = 2.
rho = -0.5
sigma = 0.5

# For restarting calculations, start the calculations as normal.
