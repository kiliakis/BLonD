import numpy as np
from input_parameters.preprocess import *
from input_parameters.general_parameters import *
import sys
from decimal import Decimal
import matplotlib.pyplot as plt
from beams.beams import *
from input_parameters.rf_parameters import *
from plots.plot_beams import *
from plots.plot_impedance import *
from plots.plot_slices import *
from plots.plot import *
from plots.plot_parameters import *
from beams.slices import *
from monitors.monitors import *
from trackers.tracker import *
import time 
from matplotlib.animation import ArtistAnimation
from beams.distributions import *
from llrf.phase_loop import *

# Beam parameters
particle_type = 'proton'
n_macroparticles = 100000
n_particles = 0


# Machine and RF parameters
radius = 25 # [m]
gamma_transition = 4.076750841  # [1]
alpha = 1 / gamma_transition**2 # [1] 
C = 2*np.pi*radius  # [m]     

n_turns = 10000

general_params = GeneralParameters(n_turns, C, alpha, 310891054.809, 
                                   particle_type)
# Cavities parameters
n_rf_systems = 1                                     
harmonic_numbers_1 = 1  # [1]  
voltage_1 = 8000  # [V]  
phi_offset_1 = 0   # [rad]
rf_params = RFSectionParameters(general_params, n_rf_systems, harmonic_numbers_1, voltage_1, phi_offset_1, omega_rf = 1.00001*2.*np.pi/general_params.t_rev[0])

my_beam = Beam(general_params, n_macroparticles, n_particles)

slices_ring = Slices(rf_params, my_beam, 200, cut_left = -0.9e-6, cut_right = 0.9e-6)

#Phase loop
configuration = {'machine': 'PSB', 'PL_gain': 0., 'RL_gain': [0.,0.], 'PL_period': 10.e-6, 'RL_period': 7}
phase_loop = PhaseLoop(general_params, rf_params, slices_ring, configuration)


#Long tracker
long_tracker = RingAndRFSection(rf_params, my_beam, periodicity = 'Off', PhaseLoop = phase_loop)

full_ring = FullRingAndRF([long_tracker])

distribution_options = {'type': 'gaussian', 'bunch_length': 200.e-9, 'density_variable': 'density_from_J'}

matched_from_distribution_density(my_beam, full_ring, distribution_options)
slices_ring.track()
long_tracker = RingAndRFSection(rf_params, my_beam, periodicity = 'Off', PhaseLoop = phase_loop)

#Monitor
bunch_monitor = BunchMonitor(general_params, rf_params, my_beam, '../output_files/TC10_output_data',
                 Slices = slices_ring, PhaseLoop = phase_loop)


#Plots
format_options = {'dirname': '../output_files/TC10_fig'}
plots = Plot(general_params, rf_params, my_beam, 1000, 10000, -0.9e-6, 0.9e-6, -1.e6, 1.e6, 
             separatrix_plot= True, Slices = slices_ring, format_options = format_options, h5file = '../output_files/TC10_output_data', PhaseLoop = phase_loop)


# Accelerator map
map_ = [long_tracker] + [slices_ring] + [bunch_monitor] + [plots] 

#phase_loop.reference += 0.00001

for i in range(1, n_turns+1):
    
    t0 = time.clock()
    for m in map_:
        m.track()   
    slices_ring.track_cuts()   
    #print time.clock()-t0
    if (i % 100 == 0): 
        print "Time step %d" %i
        print "    Radial error %.4e" %(phase_loop.drho)
        print "    Radial error, accum %.4e" %(phase_loop.drho_int)
        print "    Radial loop frequency correction %.4e 1/s" %(phase_loop.domega_RF)
        print "    RF phase %.4f rad" %(rf_params.phi_RF[0,i])
        print "    RF frequency %.6e 1/s" %(rf_params.omega_RF[0,i])
        print "    Tracker phase %.4f rad" %(long_tracker.phi_RF[0,i])
        print "    Tracker frequency %.6e 1/s" %(long_tracker.omega_RF[0,i])
    

        
print 'DONE'
