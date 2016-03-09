/*
 * general_parameters.h
 *
 *  Created on: Mar 8, 2016
 *      Author: kiliakis
 */

#ifndef INPUT_PARAMETERS_GENERALPARAMETERS_H_
#define INPUT_PARAMETERS_GENERALPARAMETERS_H_

#include "../includes/configuration.h"
#include <vector>
#include "../includes/utilities.h"
#include <math.h>
#include <numeric>
#include <string.h>
#include "../includes/constants.h"

enum particle_type {
	proton, electron, user_input, none
};

class GeneralParameters {

private:
	void eta_generation();
	void _eta0();
	void _eta1();
	void _eta2();
public:
	int n_sections;
	particle_type particle, particle_2;
	int n_turns;
	ftype mass, mass2;
	ftype charge, charge2;
	ftype cumulative_times;
	ftype *alpha;
	ftype *momentum;
	int alpha_order;
	ftype* ring_length;
	ftype ring_circumference;
	ftype ring_radius;
	ftype *beta;
	ftype *gamma;
	ftype *energy;
	ftype *kin_energy;
	ftype* cycle_time;
	ftype* f_rev, t_rev, omega_rev;
	// TODO what is the type of f_rev, t_rev, omega_rev??
	// probably they are arrays
	ftype *eta_0, *eta_1, *eta_2;

	GeneralParameters(const int n_turns, ftype* ring_length, ftype *alpha,
			const int alpha_order, ftype *momentum,
			const particle_type particle, ftype user_mass = 0,
			ftype user_charge = 0, particle_type particle2 = none,
			ftype user_mass_2 = 0, ftype user_charge_2 = 0,
			int number_of_sectrions = 1);

};

void GeneralParameters::eta_generation() {
	_eta0();
	if (alpha_order > 0)
		_eta1();
	if (alpha_order > 1)
		_eta2();
	if (alpha_order > 2)
		dprintf(
				"WARNING: Momentum compaction factor is implemented only up to 2nd order");
}

void GeneralParameters::_eta0() {
	eta_0 = new ftype[n_sections * (n_turns + 1)];
	// TODO expand this approach to all the code
	for (int i = 0; i < n_sections * (n_turns + 1); ++i) {
		int j = i / (n_turns + 1);
		eta_0[i] = alpha[j] - 1 / (gamma[i] * gamma[i]);
	}
	//dprintf("eta_0[0] = %lf\n", eta_0[0]);

}

void GeneralParameters::_eta1() {
	eta_1 = new ftype[n_sections * (n_turns + 1)];
	for (int i = 0; i < n_sections * (n_turns + 1); ++i) {
		int j = i / (n_turns + 1);
		eta_1[i] = 3 * beta[i] * beta[i] / (2 * gamma[i] * gamma[i])
				+ alpha[j + 1] - alpha[j] * eta_0[i];
	}

}

void GeneralParameters::_eta2() {
	eta_2 = new ftype[n_sections * (n_turns + 1)];
	for (int i = 0; i < n_sections * (n_turns + 1); ++i) {
		int j = i / (n_turns + 1);
		ftype betasq = beta[i] * beta[i];
		ftype gammasq = gamma[i] * gamma[i];
		eta_1[i] = -betasq * (5 * betasq - 1) / (2 * gammasq) + alpha[j + 2]
				- 2 * alpha[j] * alpha[j + 1] + alpha[j + 1] / gammasq
				+ alpha[j] * alpha[j] * eta_0[i]
				- 3 * betasq * alpha[j] / (2 * gammasq);
	}
}

GeneralParameters::GeneralParameters(const int _n_turns, ftype* _ring_length,
		ftype* _alpha, const int _alpha_order, ftype* _momentum,
		const particle_type _particle, ftype user_mass, ftype user_charge,
		particle_type _particle2, ftype user_mass_2, ftype user_charge_2,
		int number_of_sections) {

	this->particle = _particle;
	this->particle_2 = _particle2;
	this->n_sections = number_of_sections;

	if (particle == proton) {
		mass = m_p * c * c / e;
		charge = 1;
	} else if (particle == electron) {
		mass = m_e * c * c / e;
		charge = -1;
	} else if (particle == user_input) {
		mass = user_mass;
		charge = user_charge;
	} else {
		dprintf("ERROR: Particle type not recognized!");
		exit(-1);
	}

	if (particle_2 == none) {
		;
	} else if (particle_2 == proton) {
		mass2 = m_p * c * c / e;
		charge2 = 1;
	} else if (particle == electron) {
		mass2 = m_e * c * c / e;
		charge2 = -1;
	} else if (particle == user_input) {
		mass2 = user_mass_2;
		charge2 = user_charge_2;
	} else {
		dprintf("ERROR: Second particle type not recognized!");
		exit(-1);
	}

	this->n_turns = _n_turns;

	// We don't have tuple momentum == we don't need cumulative_times
	// assuming that momentum is a 2d array
	//this->momentum = new ftype[n_sections * (n_turns + 1)];

	//this->momentum = (ftype *) malloc(
	//		sizeof(ftype) * (n_turns + 1) * n_sections);

	//memcpy(this->momentum, _momentum,
	//		(n_turns + 1) * n_sections * sizeof(ftype));

	this->momentum = _momentum;

	this->alpha_order = _alpha_order;

	//this->alpha = new ftype[n_sections * (alpha_order)];
	//this->alpha = (ftype *) malloc(sizeof(ftype) * n_sections * alpha_order);

	//memcpy(this->alpha, _alpha, (alpha_order) * n_sections * sizeof(ftype));

	this->alpha = _alpha;

	this->ring_length = new ftype[n_sections];
	//this->ring_length = (ftype *) malloc(sizeof(ftype) * n_sections);
	memcpy(this->ring_length, _ring_length, sizeof(ftype) * n_sections);

	this->ring_circumference = std::accumulate(&ring_length[0],
			&ring_length[n_sections], 0);
	this->ring_radius = ring_circumference / (2 * pi);

	// TODO check what is happening with alpha and momentum as well
	// There is no way this can happen the way we have coded this
	/*
	 if ((n_sections != ring_length.size())) {
	 dprintf(
	 "ERROR: Number of sections, ring length, alpha, and/or momentum data do not match!");
	 exit(-1);
	 }
	 */
	if (n_sections > 1) {
		// TODO do some things inside here
		// Should ask danilo about this
		// Danilo told me we could skip this for now
	}

	this->gamma = new ftype[n_sections * (n_turns + 1)];
	this->beta = new ftype[n_sections * (n_turns + 1)];
	this->energy = new ftype[n_sections * (n_turns + 1)];
	this->kin_energy = new ftype[n_sections * (n_turns + 1)];

	ftype masssq = mass * mass;

	for (int i = 0; i < n_sections * (n_turns + 1); ++i) {
		ftype momentumsq = momentum[i] * momentum[i];
		this->beta[i] = sqrt(1 / (1 + (masssq / momentumsq)));
		this->gamma[i] = sqrt(1 + (momentumsq / masssq));
		this->energy[i] = sqrt(masssq + momentumsq);
		this->kin_energy[i] = energy[i] - mass;
	}

	//TODO assuming momentum is a 2d array

	if (alpha_order > 3) {
		dprintf(
				"WARNING: Momentum compaction factor is implemented only up to 2nd order");
		alpha_order = 3;
	}

	eta_generation();
}

#endif /* INPUT_PARAMETERS_GENERALPARAMETERS_H_ */
