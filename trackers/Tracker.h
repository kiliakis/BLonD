/*
 * Tracker.h
 *
 *  Created on: Mar 10, 2016
 *      Author: kiliakis
 */

#include "../input_parameters/GeneralParameters.h"
#include "../beams/Beams.h"
#include "../includes/utilities.h"
#include "sin.h"

#ifndef TRACKERS_TRACKER_H_
#define TRACKERS_TRACKER_H_

enum solver_type {
	simple, full
};

class RingAndRfSection {

private:

public:
	bool *indices_right_outside;
	bool *indices_inside_frame;
	bool *indices_left_outside;
	ftype *insiders_dE;
	ftype *left_outsiders_dE;
	ftype *left_outsiders_dt;
	void set_periodicity();
	// TODO do we actually need these parameters? (beam_dE, beam_dt)
	void kick(ftype *beam_dt, ftype * beam_dE, int index);
	void drift(ftype *beam_dt, ftype * beam_dE, int index);
	void track();
	void orizontal_cut();
	RingAndRfSection(GeneralParameters *gp, RfParameters *rf_params,
			Beams *beam, solver_type solver = simple, ftype *PhaseLoop = NULL,
			ftype * NoiseFB = NULL, bool periodicity = false, ftype *dE_max =
					NULL, bool rf_kick_interp = false, ftype* Slices = NULL,
			ftype * TotalInducedVoltage = NULL, int n_threads = 1);
	RfParameters *rf_params;
	GeneralParameters *gp;
	Beams *beam;
	solver_type solver;
	ftype *PhaseLoop;
	ftype * NoiseFB;
	bool periodicity;
	ftype *dE_max;
	bool rf_kick_interp;
	ftype* Slices;
	ftype * TotalInducedVoltage;
	int n_threads;
	ftype* acceleration_kick;
};

inline void RingAndRfSection::kick(ftype *beam_dt, ftype * beam_dE, int index) {
	// KICK
	int k = index;
	for (int j = 0; j < rf_params->n_rf; j++) {
		for (int i = 0; i < beam->n_macroparticles; i++) {
			beam_dE[i] += rf_params->voltage[k]
					* vdt::fast_sin(
							rf_params->omega_RF[k] * beam_dt[i]
									+ rf_params->phi_RF[k]);
		}
		k += beam->n_macroparticles;
	}

	// SYNCHRONOUS ENERGY CHANGE
	for (int i = 0; i < beam->n_macroparticles; i++)
		beam_dE[i] += acceleration_kick[index];

}

inline void RingAndRfSection::drift(ftype *beam_dt, ftype * beam_dE,
		int index) {

	double T = gp->t_rev[index] * rf_params->length_ratio;

	double beta = rf_params->beta(index);
	double gamma = rf_params->gamma(index);
	double energy = rf_params->energy(index);

	if (solver == simple) {
		double coeff = rf_params->eta_0(index) / (beta * beta * energy);

		for (int i = 0; i < beam->n_macroparticles; i++)
			beam_dt[i] += T * coeff * beam_dE[i];
	} else {
		const double coeff = 1. / (beta * beta * energy);
		const double eta0 = rf_params->eta_0(index) * coeff;
		const double eta1 = rf_params->eta_1(index) * coeff * coeff;
		const double eta2 = rf_params->eta_2(index) * coeff * coeff * coeff;

		if (alpha_order == 1)
			for (int i = 0; i < beam->n_macroparticles; i++)
				beam_dt[i] += T * (1. / (1. - eta0 * beam_dE[i]) - 1.);
		else if (alpha_order == 2)
			for (int i = 0; i < beam->n_macroparticles; i++)
				beam_dt[i] += T
						* (1.
								/ (1. - eta0 * beam_dE[i]
										- eta1 * beam_dE[i] * beam_dE[i]) - 1.);
		else
			for (int i = 0; i < beam->n_macroparticles; i++)
				beam_dt[i] += T
						* (1.
								/ (1. - eta0 * beam_dE[i]
										- eta1 * beam_dE[i] * beam_dE[i]
										- eta2 * beam_dE[i] * beam_dE[i]
												* beam_dE[i]) - 1.);
	}
}

inline void RingAndRfSection::track() {
	if (periodicity) {
		// Change reference of all the particles on the right of the current
		// frame; these particles skip one kick and drift
		for (int i = 0; i < beam->n_macroparticles; ++i) {
			if (indices_right_outside[i])
				beam->dt[i] -= gp->t_rev[rf_params->counter + 1];
		}
		// Synchronize the bunch with the particles that are on the right of
		// the current frame applying kick and drift to the bunch; after that
		// all the particle are in the new updated frame

		kick(beam->insiders_dt, insiders_dE, rf_params->counter);
		drift(beam->insiders_dt, insiders_dE, rf_params->counter + 1);
		// TODO why should we maintain both arrays? (insiders_dt, dt)
		for (int i = 0; i < beam->n_macroparticles; ++i) {
			if (indices_inside_frame[i]) {
				beam->dt[i] = beam->insiders_dt[i];
				beam->dE[i] = insiders_dE[i];
			}
		}
		// Check all the particles on the left of the just updated frame and
		// apply a second kick and drift to them with the previous wave after
		// having changed reference.
		// TODO left, right outsiders should also be lists
		for (int i = 0; i < beam->n_macroparticles; ++i) {
			indices_left_outside[i] = beam->dt[i] < 0;
		}

	} else {
		kick(beam->dt, beam->dE, rf_params->counter);
		drift(beam->dt, beam->dE, rf_params->counter + 1);
	}
	if (dE_max != NULL)
		orizontal_cut();
	rf_params->counter++;
	// TODO I think there is no need to make any updates, counter is enough
}

inline void RingAndRfSection::orizontal_cut() {
	// TODO resizing the array would be very expensive
	// We should think of a new data structure or way of doing this
	// Maybe the use of lists for beam-dE,dt would be a good solution
}

RingAndRfSection::RingAndRfSection(GeneralParameters *_gp,
		RfParameters *_rf_params, Beams *_beam, solver_type _solver = simple,
		ftype *_PhaseLoop = NULL, ftype * _NoiseFB = NULL, bool _periodicity =
				false, ftype *_dE_max = NULL, bool _rf_kick_interp = false,
		ftype* _Slices = NULL, ftype * _TotalInducedVoltage = NULL,
		int _n_threads = 1) {
	this->gp = _gp;
	this->rf_params = _rf_params;
	this->beam = _beam;
	this->solver = _solver;
	this->PhaseLoop = _PhaseLoop;
	this->NoiseFB = _NoiseFB;
	this->periodicity = _periodicity;
	this->dE_max = _dE_max;
	this->rf_kick_interp = _rf_kick_interp;
	this->Slices = _Slices;
	this->TotalInducedVoltage = _TotalInducedVoltage;
	this->n_threads = _n_threads;
	this->indices_left_outside = new ftype[beam->n_macroparticles];
	this->left_outsiders_dE = new ftype[beam->n_macroparticles];
	this->left_outsiders_dt = new ftype[beam->n_macroparticles];
//TODO fill unused eta arrays with zeros

	this->acceleration_kick = new ftype[rf_params->n_rf * (gp->n_turns)];
	for (int i = 0; i < rf_params->n_rf * gp->n_turns; ++i) {
		acceleration_kick[i] = -rf_params->E_increment[i];
	}

	if (solver != simple && solver != full) {
		dprintf(
				"ERROR: Choice of longitudinal solver not recognized! Aborting...");
		exit(-1);
	}

	if (gp->alpha_order > 1) {
		solver = full;
	}

	if (periodicity) {
		// TODO why not using this function anyway?
		set_periodicity();
		// TODO Here we have a bad technique, we have to deallocate t_rev and allocate it again
		// TODO we will either do this by using vectors or not do it at all if there is no actual need
		gp->t_rev.push_back(gp->t_rev.back());
	}

}

void RingAndRfSection::set_periodicity() {
	int a = 0;
	for (int i = 0; i < beam->n_macroparticles; i++)
		a += beam->dt[i] < 0;
	if (a > 0) {
		dprintf("ERROR: condition beam.dt >= 0 not true!");
		exit(-1);
	}
	this->indices_right_outside = new bool[beam->n_macroparticles];
// TODO re-thing the need of this array
	this->indices_inside_frame = new bool[beam->n_macroparticles];
// TODO assuming that beams->insiders_dt is already zeroed
// TODO assuming that insiders_dE is already zeroed
// TODO why this is not in beams as well?
	this->insiders_dE = new ftype[beam->n_macroparticles];
	for (int i = 0; i < beam->n_macroparticles; i++) {
		if (beam->dt[i] > gp->t_rev[rf_params->counter + 1]) {
			indices_right_outside[i] = true;
			indices_inside_frame[i] = false;
		} else {
			indices_inside_frame[i] = true;
			indices_inside_frame[i] = false;
		}
	}
}

#endif /* TRACKERS_TRACKER_H_ */
