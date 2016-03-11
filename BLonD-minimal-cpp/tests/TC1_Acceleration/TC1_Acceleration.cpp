/*
 * TC1_Acceleration.cpp
 *
 *  Created on: Mar 9, 2016
 *      Author: kiliakis
 */

#include "../../includes/utilities.h"
#include "../../input_parameters/GeneralParameters.h"
#include "../../input_parameters/RfParameters.h"
#include "../../beams/Beams.h"
#include "../../trackers/Tracker.h"

#include <stdio.h>

// Simulation parameters --------------------------------------------------------
// Bunch parameters
const long N_b = 1e9;           // Intensity
const long N_p = 10000;         // Macro-particles
const ftype tau_0 = 0.4e-9;          // Initial bunch length, 4 sigma [s]

// Machine and RF parameters
const ftype C = 26658.883;          // Machine circumference [m]
const long p_i = 450e9;          // Synchronous momentum [eV/c]
const ftype p_f = 460.005e9;          // Synchronous momentum, final
const long h = 35640;          // Harmonic number
const ftype V = 6e6;          // RF voltage [V]
const ftype dphi = 0;          // Phase modulation/offset
const ftype gamma_t = 55.759505;          // Transition gamma
const ftype alpha = 1.0 / gamma_t / gamma_t;    // First order mom. comp. factor
const int alpha_order = 1;
const int n_sections = 1;
// Tracking details
const int N_t = 100000;    // Number of turns to track

// Simulation setup -------------------------------------------------------------
int main(int argc, char **argv) {
	/// initializations
	timespec begin, end;
	printf("Setting up the simulation...\n\n");
	get_time(begin);
	// TODO maybe we can change the way linspace works (take an array an initialize it)
	ftype *momentum = linspace(p_i, p_f, N_t + 1);

	ftype *alpha_array = new ftype[(alpha_order + 1) * n_sections];
	std::fill_n(alpha_array, (alpha_order + 1) * n_sections, alpha);

	ftype *C_array = new ftype[n_sections];
	C_array[0] = C;

	ftype *h_array = new ftype[n_sections * (N_t + 1)];
	std::fill_n(h_array, (N_t + 1) * n_sections, h);

	ftype *V_array = new ftype[n_sections * (N_t + 1)];
	std::fill_n(V_array, (N_t + 1) * n_sections, V);

	ftype *dphi_array = new ftype[n_sections * (N_t + 1)];
	std::fill_n(dphi_array, (N_t + 1) * n_sections, dphi);

	// printf("alpha = %lf\n", alpha);

	// TODO variables must be in the correct format (arrays for all)
	GeneralParameters *general_params = new GeneralParameters(N_t, C_array,
			alpha_array, alpha_order, momentum, proton);

	// dump(general_params.gamma, N_t + 1, "gamma");
	//printf("eta_0[0] = %.8lf\n", general_params->eta_0[0]);
	//printf("eta_0[last] = %.8lf\n", general_params->eta_0[N_t]);

	// TODO maybe general_params, beam, and RfParameters could be global?

	Beams *beam = new Beams(general_params, N_p, N_b);
	RfParameters *rf_params = new RfParameters(general_params, beam, n_sections,
			h_array, V_array, dphi_array);
	//dump(rf_params->omega_RF, n_sections * (N_t + 1), "omega_RF");

	RingAndRfSection *long_tracker = new RingAndRfSection(general_params,
			rf_params, beam);
	//dump(long_tracker->acceleration_kick, N_p, "acc_kick");
	//dump(rf_params->E_increment, n_sections * (N_t), "E_increment");
	//dump(rf_params->Qs, n_sections * (N_t + 1), "Qs");
	//dump(general_params.eta_0, N_t + 1, "eta_0");

	for (int i = 1; i < N_t + 1; ++i) {
		//printf("step %d\n", i);
		//dump(beam->dE, beam->n_macroparticles, "beam->dE");
		long_tracker->track();
		beam->losses_longitudinal_cut(0, 2.5e-9);
	}
	get_time(end);
	print_time("Total Simulation Time", begin, end);

	//printf("step %d\n", N_t + 1);
	dump(beam->dE, 5, "beam->dE");
	//dump(beam->dt, beam->n_macroparticles, "beam->dt");
	printf("Done!");

}

/*
 long_tracker = RingAndRFSection(rf_params, beam)

 // Accelerator map
 map_ = [long_tracker]
 print "Map set"
 print ""

 // Tracking ---------------------------------------------------------------------
 for i in range(1, N_t+1):

 // Plot has to be done before tracking (at least for cases with separatrix)
 if (i % dt_plt) == 0:
 print "Outputting at time step %d..." %i
 print "   Beam momentum %.6e eV" %beam.momentum
 print "   Beam gamma %3.3f" %beam.gamma
 print "   Beam beta %3.3f" %beam.beta
 print "   Beam energy %.6e eV" %beam.energy
 print "   Four-times r.m.s. bunch length %.4e s" %(4.*beam.sigma_dt)
 print ""

 // Track
 for m in map_:
 m.track()

 // Define losses according to separatrix and/or longitudinal position
 beam.losses_longitudinal_cut(0., 2.5e-9)

 print "Done!"
 print ""
 */
