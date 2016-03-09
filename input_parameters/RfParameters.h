/*
 * RfParameters.h
 *
 *  Created on: Mar 9, 2016
 *      Author: kiliakis
 */

#ifndef INPUT_PARAMETERS_RFPARAMETERS_H_
#define INPUT_PARAMETERS_RFPARAMETERS_H_

#include "../input_parameters/GeneralParameters.h"
#include "../includes/utilities.h"

class RfParameters {
public:
	RfParameters(GeneralParameters *_gp, int _n_rf, ftype *_harmonic,
			ftype *_voltage, ftype *_phi_offset, ftype* _phi_noise,
			ftype * _omega_rf, int _section_index);
	ftype *E_increment;
	ftype phi_s;
	ftype *Qs;
	ftype *omega_s0;
private:
	ftype eta_0(const int i);
	ftype eta_1(const int i);
	ftype eta_2(const int i);
	int sign_eta_0(const int i);
	// TODO assume input_value is an array
	// that is why we don't have any input_check function
	int counter;
	GeneralParameters *gp;
	int n_rf;
	//int n_turns;
	ftype *harmonic;
	ftype *voltage;
	ftype *phi_offset;
	ftype *phi_noise;
	ftype *omega_rf;
	int section_index;
};

/*
 :How to use RF programs:

 - For 1 RF system and constant values of V, h or phi, just input the single value
 - For 1 RF system and varying values of V, h or phi, input an array of n_turns values
 - For several RF systems and constant values of V, h or phi, input lists of single values
 - For several RF systems and varying values of V, h or phi, input lists of arrays of n_turns values
 */
// TODO the input will have to be a single array in any case
// RfParameters == RfSectionParameters
// completely removed accelerating_systems
RfParameters::RfParameters(GeneralParameters *_gp, int _n_rf, ftype *_harmonic,
		ftype *_voltage, ftype *_phi_offset, ftype* _phi_noise = NULL,
		ftype * _omega_rf = NULL, int _section_index = 1) {
	this->counter = 0;
	this->gp = _gp;
	this->section_index = _section_index - 1;
	this->n_rf = _n_rf;
	//this->n_turns = gp->n_turns;
	this->harmonic = _harmonic;
	this->voltage = _voltage;
	this->phi_offset = _phi_offset;
	this->phi_noise = _phi_noise;
	this->omega_rf = _omega_rf;
	// TODO how to initialize this phi_s
	this->phi_s = 0;

	// wiped out all the imports
	this->E_increment = new ftype[n_rf * (gp->n_turns)];
	// Don't have n_turns +1 cause of the way np.diff works
	int k = 0;
	for (int i = 0; i < n_rf; ++i) {
		for (int j = 0; j < gp->n_turns; ++j) {
			E_increment[k] = gp->energy[k + 1] - gp->energy[k];
			k++;
		}
		k++;
	}
	// to use all these eta and sign_eta_0, there are functions that have been implemented
	// you just have to specify the number of turn you want
	// TODO no pre - processing has been done
	this->Qs = new ftype[n_rf * (gp->n_turns + 1)];
	k = 0;
	for (int i = 0; i < n_rf; ++i) {
		for (int j = 0; j < gp->n_turns + 1; ++j) {
			Qs[i] = sqrt(
					harmonic[j] * gp->charge * voltage[j]
							* abs(eta_0(k) * cos(phi_s))
							/ (2 * pi * gp->beta[k] * gp->beta[k]
									* gp->energy[k]));
			k++;
		}
	}

	this->omega_s0 = new ftype[n_rf * (gp->n_turns + 1)];
	for (int i = 0; i < n_rf * (gp->n_turns + 1); ++i) {
		this->omega_s0[i] = Qs[i] * gp->omega_rev;
	}

	// TODO this where I should continue
}

/*
 

 #: *Design RF frequency of the RF systems in the station [Hz]*        
 self.omega_RF_d = 2.*np.pi*self.beta*c*self.harmonic/ \
                          (self.ring_circumference)
 
 #: *Initial, actual RF frequency of the RF systems in the station [Hz]*
 if omega_rf == None:
 self.omega_RF = np.array(self.omega_RF_d)                  

 #: *Initial, actual RF phase of each harmonic system*
 self.phi_RF = np.array(self.phi_offset) 
 
 #: *Accumulated RF phase error of each harmonic system*
 self.dphi_RF = np.zeros(self.n_rf)
 
 #: *Accumulated RF phase error of each harmonic system*
 self.dphi_RF_steering = np.zeros(self.n_rf)
 
 self.t_RF = 2*np.pi / self.omega_RF[0]
 
 
 
 
 def eta_tracking(self, beam, counter, dE):
 '''
 *The slippage factor is calculated as a function of the energy offset
 (dE) of the beam particle. By definition, the slippage factor in ith 
 order is:*
 
 .. math:: 
 \\eta(\\delta) = \\sum_{i}(\\eta_i \\, \\delta^i) = \\sum_{i} \\left(\\eta_i \\, \\left[ \\frac{\\Delta E}{\\beta_s^2 E_s} \\right]^i \\right)
 
 '''
 
 if self.alpha_order == 1:
 return self.eta_0[counter]
 else:
 eta = 0
 delta = dE/(beam.beta**2 * beam.energy)
 for i in xrange( self.alpha_order ):
 eta_i = getattr(self, 'eta_' + str(i))[counter]
 eta  += eta_i * (delta**i)
 return eta  

 */

ftype RfParameters::eta_0(const int i) {
	return gp->eta_0[section_index * (gp->n_turns + 1) + i];
}

ftype RfParameters::eta_1(const int i) {
	return gp->eta_1[section_index * (gp->n_turns + 1) + i];
}

ftype RfParameters::eta_2(const int i) {
	return gp->eta_2[section_index * (gp->n_turns + 1) + i];
}

int RfParameters::sign_eta_0(const int i) {
	if (gp->eta_0[section_index * (gp->n_turns + 1) + i] > 0)
		return 1;
	else if (gp->eta_0[section_index * (gp->n_turns + 1) + i] == 0)
		return 0;
	else
		return -1;
}

#endif /* INPUT_PARAMETERS_RFPARAMETERS_H_ */
