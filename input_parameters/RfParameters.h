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
	ftype *omega_RF_d;
	ftype *phi_rf;
	ftype *dphi_RF;
	ftype *dphi_RF_steering;
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
	frype *omega_RF;
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

	this->omega_RF_d = new ftype[n_rf*(gp->n_turns+1)];
	for (int i = 0; i < n_rf * (gp->n_turns + 1); ++i) {
		this->omega_RF_d[i] = 2*pi*gp->beta[i] * harmonic[i]/ gp->ring_circumference;
	}

	this->omega_RF = new ftype[n_rf *(gp->n_turns +1)];
	if (omega_rf == NULL){
		std::copy(std::begin(omega_RF_d), std::end(omega_RF_d), omega_RF);
	}

	this->phi_RF = new ftype[n_rf *(gp->n_turns +1)];
	std::copy(std::begin(phi_offset), std::end(phi_offset), phi_RF);

	this->dphi_RF = new ftype[n_rf];
	std::fill_n(dphi_RF, n_rf, 0);

	this->dphi_RF_steering = new ftype[n_rf];
	std::fill_n(dphi_RF_steering, n_rf, 0);

	// TODO probably this is an array
	this->t_RF = 2*pi / omega_RF[0];
}
	
	// TODO what is this beam.beta, beam.energy thing?
	// TODO maybe I should add some functions about this
ftype RfParameters::eta_tracking(Beams beam, int counter,ftype dE){
	if (gp->alpha_order == 1)
		return eta_0(counter);
	else{
		ftype eta = 0;
		ftype delta = dE / ((beam->gp->beta[0][0])*(beam->gp->beta[0][0])*
					beam->gp->energy[0][0]);
	eta += eta_0(counter) * 1;
	if (alpha_order > 0)
		eta += eta_1(counter) * delta;
	if (alpha_order > 1)
		eta += eta_2(counter) * delta*delta;
	if (alpha_order > 2)
		dprintf(
				"WARNING: Momentum compaction factor is implemented only up to 2nd order");
	}
	return eta;

}

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
