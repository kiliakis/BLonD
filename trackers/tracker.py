# Copyright 2015 CERN. This software is distributed under the
# terms of the GNU General Public Licence version 3 (GPL Version 3), 
# copied verbatim in the file LICENCE.md.
# In applying this licence, CERN does not waive the privileges and immunities 
# granted to it by virtue of its status as an Intergovernmental Organization or
# submit itself to any jurisdiction.
# Project website: http://blond.web.cern.ch/

'''
**Module containing all the elements to track the RF frequency and phase and the
beam in phase space.**

:Authors:  **Helga Timko**, **Alexandre Lasheen**, **Danilo Quartullo**
'''

from __future__ import division
import numpy as np
from scipy.integrate import cumtrapz
import ctypes
from setup_cpp import libfib
from scipy.constants import c
import matplotlib.pyplot as plt




class RingAndRFSection(object):
    '''
    *Definition of an RF station and part of the ring until the next station, 
    see figure.*
    
    .. image:: ring_and_RFstation.png
        :align: center
        :width: 600
        :height: 600
        
    *The time step is fixed to be one turn, but the tracking can consist of 
    multiple RingAndRFSection objects. In this case, the user should make sure 
    that the lengths of the stations sum up exactly to the circumference or use
    the FullRingAndRF object in order to let the code pre-process the 
    parameters. Each RF station may contain several RF harmonic systems which 
    are considered to be in the same location. First, the energy kick of the RF
    station is applied, and then the particle arrival time to the next station
    is updated. The change in RF phase and frequency due to control loops is 
    tracked as well.*
    '''
        
    def __init__(self, RFSectionParameters, Beam, solver = 'simple', 
                 PhaseLoop = None, NoiseFB = None, periodicity = False, dE_max = None, rf_kick_interp=False, Slices=None, TotalInducedVoltage=None, n_threads=1):
        
        #: *Import of RFSectionParameters object*
        self.rf_params = RFSectionParameters

        #: *Import of Beam object*
        self.beam = Beam
        
        ### Import RF section parameters #######################################
        #: *Import section index (from RFSectionParameters)*        
        self.section_index = RFSectionParameters.section_index
        
        #: *Import counter (from RFSectionParameters)*        
        self.counter = RFSectionParameters.counter 
              
        #: *Import length ratio (from RFSectionParameters)*
        self.length_ratio = RFSectionParameters.length_ratio
        
        #: *Import section length (from RFSectionParameters)* # needed for FullRingAndRF
        self.section_length = RFSectionParameters.section_length
        
        #: *Import revolution period (from GeneralParameters)*       
        self.t_rev = RFSectionParameters.t_rev

        #: *Import the number of RF systems (from RFSectionParameters)*
        self.n_rf = RFSectionParameters.n_rf
        
        #: *Import beta (from RFSectionParameters)* # needed for FullRingAndRF
        self.beta = RFSectionParameters.beta
        
        #: *Import particle charge (from RFSectionParameters)* 
        self.charge = RFSectionParameters.charge
        
        #: *Import RF harmonic number program (from RFSectionParameters)*
        self.harmonic = RFSectionParameters.harmonic 
               
        #: *Import RF voltage program [V] (from RFSectionParameters)*
        self.voltage = RFSectionParameters.voltage  
           
        #: *Import RF phase noise [rad] (from RFSectionParameters)*
        self.phi_noise = RFSectionParameters.phi_noise
        
        #: *Import RF phase [rad] (from RFSectionParameters)*
        self.phi_RF = RFSectionParameters.phi_RF
        
        #: *Import phi_s [rad] (from RFSectionParameters)* # needed for FullRingAndRF
        self.phi_s = RFSectionParameters.phi_s
        
        #: *Import actual RF frequency [1/s] (from RFSectionParameters)*
        self.omega_RF = RFSectionParameters.omega_RF
        
        #: *Slippage factor (0th order) for the given RF section*
        self.eta_0 = RFSectionParameters.eta_0
        
        #: *Slippage factor (1st order) for the given RF section*
        self.eta_1 = RFSectionParameters.eta_1
        
        #: *Slippage factor (2nd order) for the given RF section*
        self.eta_2 = RFSectionParameters.eta_2
        
        #: *Slippage factor (2nd order) for the given RF section*
        self.sign_eta_0 = RFSectionParameters.sign_eta_0
        
        #: *Import alpha order (from RFSectionParameters)*                
        self.alpha_order = RFSectionParameters.alpha_order
        
        #: *Fill unused eta arrays with zeros*
        for i in xrange( self.alpha_order, 3 ):
            setattr(self, "eta_%s" %i, np.zeros(RFSectionParameters.n_turns+1))       
        ### End of import of RF section parameters #############################
            
        #: *Synchronous energy change* :math:`: \quad - \delta E_s`
        self.acceleration_kick = - RFSectionParameters.E_increment  
        
        #: | *Choice of drift solver options*
        self.solver = solver
        if self.solver != 'simple' and self.solver != 'full':
            raise RuntimeError("ERROR: Choice of longitudinal solver not recognized! Aborting...")
            
        #: | *Set to 'full' if higher orders of eta are used*
        if self.alpha_order > 1:
            self.solver = 'full'
        
        # Set the horizontal cut
        self.dE_max = dE_max
        
        # Periodicity setting up
        self.periodicity = periodicity
        if periodicity:
            # Check the periodicity loop invariant dt>=0.
            if len(np.where(self.beam.dt<0)[0])>0:
                raise RuntimeError('ERROR: condition beam.dt >= 0 not true!')
            # Distinguish the particle inside the frame from the particles on the
            # right of the frame.
            self.indices_right_outside = np.where(self.beam.dt > self.t_rev[self.counter[0]+1])[0]
            self.indices_inside_frame = np.where(self.beam.dt < self.t_rev[self.counter[0]+1])[0]
            self.beam.insiders_dt = np.ascontiguousarray(self.beam.dt[self.indices_inside_frame])
            self.insiders_dE = np.ascontiguousarray(self.beam.dE[self.indices_inside_frame])
            self.t_rev = np.append(self.t_rev, self.t_rev[-1])
            
        # Use interpolate to apply kick
        self.n_threads = n_threads
            
    def set_periodicity(self):
        
        # Check the periodicity loop invariant dt>=0.
        if len(np.where(self.beam.dt<0)[0])>0:
            raise RuntimeError('ERROR: condition beam.dt >= 0 not true!')
        # Distinguish the particle inside the frame from the particles on the
        # right of the frame.
        self.indices_right_outside = np.where(self.beam.dt > self.t_rev[self.counter[0]+1])[0]
        self.indices_inside_frame = np.where(self.beam.dt < self.t_rev[self.counter[0]+1])[0]
        self.beam.insiders_dt = np.ascontiguousarray(self.beam.dt[self.indices_inside_frame])
        self.insiders_dE = np.ascontiguousarray(self.beam.dE[self.indices_inside_frame])        
    
    
    def kick(self, beam_dt, beam_dE, index):
        '''
        *Update of the particle energy due to the RF kick in a given RF station. 
        The kicks are summed over the different harmonic RF systems in the 
        station. The cavity phase can be shifted by the user via phi_offset.
        The main RF (harmonic[0]) has by definition phase=0 at time=0. The 
        phases of all other RF systems are defined w.r.t. to the main RF.
        The increment in energy is given by the discrete equation of motion:*
        
        .. math::
            \Delta E^{n+1} = \Delta E^n + \sum_{k=0}^{n_{\mathsf{rf}}-1}{e V_k^n \\sin{\\left(\omega_{\mathsf{rf,k}}^n \\Delta t^n + \phi_{\mathsf{rf,k}}^n \\right)}} - (E_s^{n+1} - E_s^n) 
            
        '''
        
        voltage_kick = np.ascontiguousarray(self.charge*
                                      self.voltage[:, index])
        omegaRF_kick = np.ascontiguousarray(self.omega_RF[:, index])
        phiRF_kick = np.ascontiguousarray(self.phi_RF[:, index])
        
        libfib.kick(beam_dt.ctypes.data_as(ctypes.c_void_p), 
            beam_dE.ctypes.data_as(ctypes.c_void_p), 
            ctypes.c_int(self.n_rf), voltage_kick.ctypes.data_as(ctypes.c_void_p), 
            omegaRF_kick.ctypes.data_as(ctypes.c_void_p), 
            phiRF_kick.ctypes.data_as(ctypes.c_void_p),
            ctypes.c_int(len(beam_dt)), 
            ctypes.c_double(self.acceleration_kick[index]))
        
   
    def drift(self, beam_dt, beam_dE, index):
        '''
        *Update of particle arrival time to the RF station. If only the zeroth 
        order slippage factor is given, 'simple' and 'full' solvers are 
        available. The 'simple' solver is somewhat faster. Otherwise, the solver
        is automatically 'full' and calculates the frequency slippage up to 
        second order.*
        
        *The corresponding equations are:*
        
        .. math::
            \\Delta t^{n+1} = \\Delta t^{n} + \\frac{L}{C} T_0^{n+1} \\left(\\frac{1}{1 - \\eta(\\delta^{n+1})\\delta^{n+1}} - 1\\right) \quad \\text{(full)}
            
        .. math::
            \\Delta t^{n+1} = \\Delta t^{n} + \\frac{L}{C} T_0^{n+1}\\eta_0\\delta^{n+1} \quad \\text{(simple)}
        
        '''
        
        libfib.drift(beam_dt.ctypes.data_as(ctypes.c_void_p), 
            beam_dE.ctypes.data_as(ctypes.c_void_p), 
            ctypes.c_char_p(self.solver),
            ctypes.c_double(self.t_rev[index]),
            ctypes.c_double(self.length_ratio), 
            ctypes.c_double(self.alpha_order), 
            ctypes.c_double(self.eta_0[index]), 
            ctypes.c_double(self.eta_1[index]),
            ctypes.c_double(self.eta_2[index]), 
            ctypes.c_double(self.rf_params.beta[index]), 
            ctypes.c_double(self.rf_params.energy[index]), 
            ctypes.c_int(len(beam_dt)),
            ctypes.c_int(self.n_threads))


    def track(self):
        '''
        *Tracking method for the section. Applies first the kick, then the 
        drift. Calls also RF feedbacks if applicable. Updates the counter of the
        corresponding RFSectionParameters class and the energy-related 
        variables of the Beam class.*
        '''
        
         if self.periodicity:
            
            # Change reference of all the particles on the right of the current
            # frame; these particles skip one kick and drift
            if len(self.indices_right_outside)>0:
                self.beam.dt[self.indices_right_outside] -= self.t_rev[self.counter[0]+1]
            
            # Syncronize the bunch with the particles that are on the right of
            # the current frame applying kick and drift to the bunch; after that 
            # all the particle are in the new updated frame
            self.kick(self.beam.insiders_dt, self.insiders_dE, self.counter[0])
            self.drift(self.beam.insiders_dt, self.insiders_dE, self.counter[0]+1)
            self.beam.dt[self.indices_inside_frame] = self.beam.insiders_dt
            self.beam.dE[self.indices_inside_frame] = self.insiders_dE
            
            # Check all the particles on the left of the just updated frame and 
            # apply a second kick and drift to them with the previous wave after
            # having changed reference.
            self.indices_left_outside = np.where(self.beam.dt < 0)[0]
            if len(self.indices_left_outside)>0:
                left_outsiders_dt = np.ascontiguousarray(self.beam.dt[self.indices_left_outside])
                left_outsiders_dE = np.ascontiguousarray(self.beam.dE[self.indices_left_outside])
                left_outsiders_dt += self.t_rev[self.counter[0]+1]
                self.kick(left_outsiders_dt, left_outsiders_dE, self.counter[0])
                self.drift(left_outsiders_dt, left_outsiders_dE, self.counter[0]+1)
                self.beam.dt[self.indices_left_outside] = left_outsiders_dt
                self.beam.dE[self.indices_left_outside] = left_outsiders_dE
            
            # Distinguish the particles inside the frame from the particles on the
            # right of the frame.
            self.indices_right_outside = np.where(self.beam.dt > self.t_rev[self.counter[0]+2])[0]
            self.indices_inside_frame = np.where(self.beam.dt < self.t_rev[self.counter[0]+2])[0]
            self.beam.insiders_dt = np.ascontiguousarray(self.beam.dt[self.indices_inside_frame])
            self.insiders_dE = np.ascontiguousarray(self.beam.dE[self.indices_inside_frame])
            
            # Orizzontal cut: this method really eliminates particles from the
            # code
            if self.dE_max!=None:
                itemindex = np.where(self.beam.dE > -self.dE_max)[0]
                self.beam.dt = np.ascontiguousarray(self.beam.dt[itemindex])
                self.beam.dE = np.ascontiguousarray(self.beam.dE[itemindex])
                self.beam.n_macroparticles = len(self.beam.dt)
                
        else:
            
            self.kick(self.beam.dt, self.beam.dE, self.counter[0])
            
            self.drift(self.beam.dt, self.beam.dE, self.counter[0]+1)
            
            # Orizzontal cut: this method really eliminates particles from the
            # code
            if self.dE_max!=None:
                itemindex = np.where((self.beam.dE > -self.dE_max)&(self.beam.dE < self.dE_max))[0]
                self.beam.dt = np.ascontiguousarray(self.beam.dt[itemindex])
                self.beam.dE = np.ascontiguousarray(self.beam.dE[itemindex])
                self.beam.n_macroparticles = len(self.beam.dt)
    
        # Increment by one the turn counter
        self.counter[0] += 1
        
        # Updating the beam synchronous momentum etc.
        self.beam.beta = self.rf_params.beta[self.counter[0]]
        self.beam.gamma = self.rf_params.gamma[self.counter[0]]
        self.beam.energy = self.rf_params.energy[self.counter[0]]
        self.beam.momentum = self.rf_params.momentum[self.counter[0]]
