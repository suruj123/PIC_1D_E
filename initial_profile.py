import numpy as np
import matplotlib.pyplot as plt
import scipy.sparse as sp
from scipy.sparse.linalg import spsolve	




# Simulation parameters
N         = 40000   # Number of particles
Nx        = 400     # Number of mesh cells
t         = 0       # current time of the simulation
tEnd      = 50      # time at which simulation ends
dt        = 1       # timestep
boxsize   = 50      # periodic domain [0,boxsize]
n0        = 1       # electron number density
vb        = 3       # beam velocity
vth       = 1       # beam width
A         = 0.1     # perturbation
plotRealTime = True # switch on for plotting as the simulation goes along
	
# Generate Initial Conditions
np.random.seed(42)            # set the random number generator seed
# construct 2 opposite-moving Guassian beams
pos  = np.random.rand(N,1) * boxsize  
pos = pos - boxsize*0.5
vel  = vth * np.random.randn(N,1) + vb
Nh = int(N/2)
vel[Nh:] *= -1
# add perturbation
vel *= (1 + A*np.sin(2*np.pi*pos/boxsize))
c = np.savetxt("STATE",np.c_[ pos, vel], delimiter = '\t', fmt = '%0.16f')

	
