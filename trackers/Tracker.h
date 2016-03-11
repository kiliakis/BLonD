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

// TODO
// !!!!WARNING!!!!
// we only use beam->dt, dE
// bool arrays are used every time to update only the right values!!

class RingAndRfSection {

private:

public:
	bool *indices_right_outside;
	bool *indices_inside_frame;
	bool *indices_left_outside;
	void set_periodicity();
	void kick(bool *update, int index);
	void drift(bool *update, int index);
	void track();
	void orizontal_cut();
	RingAndRfSection(GeneralParameters *gp, RfParameters *rf_params,
			Beams *beam, solver_type solver = simple, ftype *PhaseLoop = NULL,
			ftype * NoiseFB = NULL, bool periodicity = false, ftype dE_max = 0,
			bool rf_kick_interp = false, ftype* Slices = NULL,
			ftype * TotalInducedVoltage = NULL, int n_threads = 1);
	RfParameters *rf_params;
	GeneralParameters *gp;
	Beams *beam;
	solver_type solver;
	ftype *PhaseLoop;
	ftype * NoiseFB;
	bool periodicity;
	ftype dE_max;
	bool rf_kick_interp;
	ftype* Slices;
	ftype * TotalInducedVoltage;
	int n_threads;
	ftype* acceleration_kick;
};

inline void RingAndRfSection::kick(bool *update, int index) {
	// KICK
	int k = index;
	for (int j = 0; j < rf_params->n_rf; j++) {
		for (int i = 0; i < beam->n_macroparticles; i++) {
			beam->dE[i] += update[i] * rf_params->voltage[k]
					* vdt::fast_sin(
							rf_params->omega_RF[k] * beam->dt[i]
									+ rf_params->phi_RF[k]);
		}
		k += beam->n_macroparticles;
	}

	// SYNCHRONOUS ENERGY CHANGE
	for (int i = 0; i < beam->n_macroparticles; i++)
		beam->dE[i] += update[i] * acceleration_kick[index];

}

inline void RingAndRfSection::drift(bool *update, int index) {

	double T = gp->t_rev[index] * rf_params->length_ratio;

	double beta = rf_params->beta(index);
	double energy = rf_params->energy(index);

	if (solver == simple) {
		double coeff = rf_params->eta_0(index) / (beta * beta * energy);

		for (int i = 0; i < beam->n_macroparticles; i++)
			beam->dt[i] += update[i] * T * coeff * beam->dE[i];
	} else {
		const double coeff = 1. / (beta * beta * energy);
		const double eta0 = rf_params->eta_0(index) * coeff;
		const double eta1 = rf_params->eta_1(index) * coeff * coeff;
		const double eta2 = rf_params->eta_2(index) * coeff * coeff * coeff;

		if (gp->alpha_order == 1)
			for (int i = 0; i < beam->n_macroparticles; i++)
				beam->dt[i] += update[i] * T
						* (1. / (1. - eta0 * beam->dE[i]) - 1.);
		else if (gp->alpha_order == 2)
			for (int i = 0; i < beam->n_macroparticles; i++)
				beam->dt[i] +=
						update[i] * T
								* (1.
										/ (1. - eta0 * beam->dE[i]
												- eta1 * beam->dE[i]
														* beam->dE[i]) - 1.);
		else
			for (int i = 0; i < beam->n_macroparticles; i++)
				beam->dt[i] += update[i] * T
						* (1.
								/ (1. - eta0 * beam->dE[i]
										- eta1 * beam->dE[i] * beam->dE[i]
										- eta2 * beam->dE[i] * beam->dE[i]
												* beam->dE[i]) - 1.);
	}
}

inline void RingAndRfSection::track() {
	if (periodicity) {
		// Change reference of all the particles on the right of the current
		// frame; these particles skip one kick and drift
		for (int i = 0; i < beam->n_macroparticles; ++i) {
			beam->dt[i] -= indices_right_outside[i]
					* gp->t_rev[rf_params->counter + 1];
		}
		// Synchronize the bunch with the particles that are on the right of
		// the current frame applying kick and drift to the bunch; after that
		// all the particle are in the new updated frame

		kick(indices_inside_frame, rf_params->counter);
		drift(indices_inside_frame, rf_params->counter + 1);

		// find left outside particles and kick, drift them one more time
		int a = 0;
		for (int i = 0; i < beam->n_macroparticles; ++i) {
			if (beam->dt[i] < 0) {
				indices_left_outside[i] = beam->id[i] > 0;
				a++;
			} else {
				indices_left_outside[i] = false;
			}
		}
		if (a > 0) {
			// This will update only the indices_left_outside values
			//  need to test this
			for (int i = 0; i < beam->n_macroparticles; ++i) {
				beam->dt[i] += gp->t_rev[rf_params->counter + 1]
						* indices_left_outside[i];
			}
			kick(indices_left_outside, rf_params->counter);
			drift(indices_left_outside, rf_params->counter + 1);

		}
		// update inside, right outside particles

		set_periodicity();

	} else {
		kick(indices_inside_frame, rf_params->counter);
		drift(indices_inside_frame, rf_params->counter + 1);
	}
	// cut particles by zeroing their id
	// this way they will not be considered again in an update
	// TODO maybe using lists and actually removing them can be
	// faster, I need to test it (although I don't believe it)
	if (dE_max > 0)
		orizontal_cut();
	rf_params->counter++;
// TODO I think there is no need to make any updates, counter is enough
}

inline void RingAndRfSection::orizontal_cut() {
	// TODO resizing the array would be very expensive
	// We should think of a new data structure or way of doing this
	// Maybe the use of lists for beam-dE,dt would be a good solution
	// In order to cut a particle we will 0 its id
	for (int i = 0; i < beam->n_macroparticles; ++i) {
		if (beam->dE[i] < -dE_max || beam->dE[i] > dE_max)
			beam->id[i] = 0;
	}
}

RingAndRfSection::RingAndRfSection(GeneralParameters *_gp,
		RfParameters *_rf_params, Beams *_beam, solver_type _solver,
		ftype *_PhaseLoop, ftype * _NoiseFB, bool _periodicity, ftype _dE_max,
		bool _rf_kick_interp, ftype* _Slices, ftype * _TotalInducedVoltage,
		int _n_threads) {
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
	this->indices_left_outside = new bool[beam->n_macroparticles];
	this->indices_right_outside = new bool[beam->n_macroparticles];
	this->indices_inside_frame = new bool[beam->n_macroparticles];
	for (int i = 0; i < beam->n_macroparticles; ++i) {
		indices_inside_frame[i] = true;
		indices_left_outside[i] = indices_right_outside[i] = false;
	}

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
		int a = 0;
		for (int i = 0; i < beam->n_macroparticles; i++)
			a += beam->dt[i] < 0;
		if (a > 0) {
			dprintf("ERROR: condition beam.dt >= 0 not true!");
			exit(-1);
		}

		set_periodicity();
		// TODO Here we have a bad technique, we have to deallocate t_rev and allocate it again
		// TODO we will either do this by using vectors or not do it at all if there is no actual need
		gp->t_rev.push_back(gp->t_rev.back());
	}

}

inline void RingAndRfSection::set_periodicity() {

	for (int i = 0; i < beam->n_macroparticles; i++) {
		if (beam->dt[i] > gp->t_rev[rf_params->counter + 1]) {
			indices_right_outside[i] = beam->id[i] > 0;
			indices_inside_frame[i] = false;
		} else {
			indices_inside_frame[i] = beam->id[i] > 0;
			indices_inside_frame[i] = false;
		}
	}
}

#endif /* TRACKERS_TRACKER_H_ */
