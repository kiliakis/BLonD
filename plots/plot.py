
# Copyright 2015 CERN. This software is distributed under the
# terms of the GNU General Public Licence version 3 (GPL Version 3), 
# copied verbatim in the file LICENCE.md.
# In applying this licence, CERN does not waive the privileges and immunities 
# granted to it by virtue of its status as an Intergovernmental Organization or
# submit itself to any jurisdiction.
# Project website: http://blond.web.cern.ch/

'''
**Class to choose plots and customize plot layout**

:Authors: **Helga Timko**
'''

import os
import subprocess
import sys
import warnings
import matplotlib.pyplot as plt
import h5py as hp
from plots.plot_beams import *
from plots.plot_slices import *
from plots.plot_llrf import *



def fig_folder(dirname):
    '''
    Create folder where plots will be stored.
    '''
    
    # Try to create directory
    try:
        os.makedirs(dirname)
    # Check whether already exists/creation failed
    except OSError:
        if os.path.exists(dirname):
            pass
        else:
            raise



class Plot(object):
    
    def __init__(self, GeneralParameters, RFSectionParameters, Beam, dt_plot,
                 dt_bckp, xmin, xmax, ymin, ymax, xunit = 's', sampling = 1, 
                 separatrix_plot = False, histograms_plot = True, 
                 Slices = None, h5file = None, output_frequency = 1, 
                 PhaseLoop = None, LHCNoiseFB = None, format_options = None):
        '''
        Define what plots should be plotted during the simulation. Passing only
        basic objects, only phase space plot will be produced. Passing optional
        objects, plots related to those objects will be produced as well.
        For plots at a certain turn: use 'dt_plot' to set the plotting frequency 
        in units of time steps. 
        For plots as a function of time: use 'dt_bckp' to set plotting frequency
        in units of time steps.
        '''

        #: | *Import GeneralParameters*
        self.general_params = GeneralParameters
        
        #: | *Import RFSectionParameters*       
        self.rf_params = RFSectionParameters
        
        #: | *Import actual time step RFSectionParameters*
        self.tstep = RFSectionParameters.counter
        
        #: | *Import Beam*
        self.beam = Beam
        
        #: | *Plotting frequency in units of time steps*
        self.dt_plot = dt_plot
        self.dt_bckp = dt_bckp
        
        #: | *Plot limit (where applicable): minimum on x-axis [xunit]*
        self.xmin = xmin

        #: | *Plot limit (where applicable): maximum on x-axis [xunit]*
        self.xmax = xmax

        #: | *Choice of x-axis unit (where applicable): 's' or 'rad'*        
        self.xunit = xunit

        #: | *Plot limit (where applicable): minimum on y-axis [eV]*
        self.ymin = ymin

        #: | *Plot limit (where applicable): maximum on y-axis [eV]*
        self.ymax = ymax

        #: | *Sampling of large arrays (where applicable)*
        self.sampling = sampling

        #: | *Separatrix in phase space plot: 'True' or 'False'*
        self.separatix = separatrix_plot

        #: | *Histogram in phase space plot: 'True' or 'False'*
        self.histogram = histograms_plot

        #: | *Optional import of Slices*
        self.slices = Slices

        #: | *Optional import of Monitor file*
        self.h5file = h5file

        #: | *Optional sampling of monitored quantities*
        self.dt_mon = output_frequency
        
        #: | *Optional import of PhaseLoop*
        self.PL = PhaseLoop

        #: | *Optional import of LHCNoiseFB*
        self.noiseFB = LHCNoiseFB
        
        # Set plotting format
        self.set_format(format_options)
        
        # Track at initialisation
        self.track()          


    
    def set_format(self, format_options):
        '''
        Initialize plot folder and custom plot formatting. For more options, see
        
        http://matplotlib.org/1.3.1/users/customizing.html
        '''

        if not format_options.has_key('dirname'):  
            self.dirname = 'fig'
        else: 
            self.dirname = format_options['dirname'] 
            
        if not format_options.has_key('linewidth'):  
            self.lwidth = 2
        else: 
            self.lwidth = format_options['linewidth'] 

        if not format_options.has_key('linestyle'):  
            self.lstyle = '-'
        else: 
            self.lstyle = format_options['linestyle'] 

        if not format_options.has_key('markersize'):  
            self.msize = 6
        else: 
            self.msize = format_options['markersize'] 

        if not format_options.has_key('alpha'):  
            self.alpha = 0.05
        else: 
            self.alpha = format_options['alpha'] 

        if not format_options.has_key('labelsize'):  
            self.lsize = 18
        else: 
            self.lsize = format_options['labelsize'] 

        if not format_options.has_key('fontfamily'):  
            self.ffamily = 'sans-serif'
        else: 
            self.ffamily = format_options['fontfamily'] 
            
        if not format_options.has_key('fontweight'):  
            self.fweight = 'normal'
        else: 
            self.fweight = format_options['fontweight'] 

        if not format_options.has_key('dpi'):  
            self.dpi = 100
        else: 
            self.dpi = format_options['dpi'] 
        
        # Directory where longitudinal_plots will be stored
        fig_folder(self.dirname)
        
        # Ticksize
        self.tsize = self.lsize - 2
        
        # Set size of x- and y-grid numbers
        plt.rc('xtick', labelsize=self.tsize) 
        plt.rc('ytick', labelsize=self.tsize)
        
        # Set x- and y-grid labelsize and weight
        plt.rc('axes', labelsize=self.lsize)
        plt.rc('axes', labelweight=self.fweight)

        # Set linewidth for continuous, markersize for discrete plotting
        plt.rc('lines', linewidth=self.lwidth, markersize=self.msize)
        
        # Set figure resolution, font
        plt.rc('figure', dpi=self.dpi)  
        plt.rc('savefig', dpi=self.dpi)  
        plt.rc('font', family=self.ffamily)  

        
    def track(self):
        '''
        Plot in certain time steps and depending on imported objects
        '''
        
        # Snapshot-type plots
        if (self.tstep[0] % self.dt_plot) == 0:
            
            plot_long_phase_space(self.general_params, self.rf_params, 
                                  self.beam, self.xmin, self.xmax, self.ymin, 
                                  self.ymax, self.xunit, 
                                  sampling = self.sampling, 
                                  separatrix_plot = self.separatix, 
                                  histograms_plot = self.histogram, 
                                  dirname = self.dirname, alpha = self.alpha)
            
            if self.slices:
                plot_beam_profile(self.slices, self.tstep[0], 
                                  style = self.lstyle, dirname = self.dirname)

                self.slices.beam_spectrum_generation(self.slices.n_slices)
                plot_beam_spectrum(self.slices, self.tstep[0], 
                                   style = self.lstyle, dirname = self.dirname)
        
        # Plots as a function of time        
        if (self.tstep[0] % self.dt_bckp) == 0 and self.h5file:
            
            h5data = h5py.File(self.h5file + '.h5', 'r')
            plot_bunch_length_evol(self.rf_params, h5data, 
                                   output_freq = self.dt_mon, 
                                   dirname = self.dirname)
            
            if self.slices and self.slices.fit_option == 'gaussian':
                    plot_bunch_length_evol_gaussian(self.rf_params, self.slices,
                                                    h5data, 
                                                    output_freq = self.dt_mon, 
                                                    dirname = self.dirname)
            plot_position_evol(self.rf_params, h5data, 
                               output_freq = self.dt_mon, 
                               style = self.lstyle, dirname = self.dirname)
            plot_energy_evol(self.rf_params, h5data, output_freq = self.dt_mon, 
                             style = self.lstyle, dirname = self.dirname)
            plot_COM_motion(self.general_params, self.rf_params, h5data,
                            output_freq = self.dt_mon, dirname = self.dirname)
            plot_transmitted_particles(self.rf_params, h5data, 
                                       output_freq = self.dt_mon, 
                                       style = self.lstyle, 
                                       dirname = self.dirname)
                    
            if self.PL:
                plot_PL_RF_freq(self.rf_params, h5data, 
                                output_freq = self.dt_mon,
                                dirname = self.dirname)
                plot_PL_RF_phase(self.rf_params, h5data, 
                                 output_freq = self.dt_mon,
                                 dirname = self.dirname)
                plot_PL_bunch_phase(self.rf_params, h5data, 
                                    output_freq = self.dt_mon, 
                                    dirname = self.dirname)
                plot_PL_phase_corr(self.rf_params, h5data, 
                                   output_freq = self.dt_mon, 
                                   dirname = self.dirname)
                plot_PL_freq_corr(self.rf_params, h5data, 
                                  output_freq = self.dt_mon, 
                                  dirname = self.dirname)
                plot_RF_phase_error(self.rf_params, h5data, 
                                    output_freq = self.dt_mon, 
                                    dirname = self.dirname)
                plot_RL_radial_error(self.rf_params, h5data, 
                                     output_freq = self.dt_mon, 
                                     dirname = self.dirname)
            
            if self.noiseFB:
                plot_LHCNoiseFB(self.rf_params, self.noiseFB, h5data, 
                                output_freq = self.dt_mon, 
                                dirname = self.dirname)
                plot_LHCNoiseFB_FWHM(self.rf_params, self.noiseFB, h5data, 
                                     output_freq = self.dt_mon, 
                                     dirname = self.dirname)

            h5data.close()


    def reset_frame(self, xmin, xmax, ymin, ymax):
        
        self.xmin = xmin
        self.xmax = xmax
        self.ymin = ymin
        self.ymax = ymax



