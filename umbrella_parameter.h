#ifndef UMBRELLA_PARAMETER_H
#define UMBRELLA_PARAMETER_H

#include <mpi.h>
#include <lammps.h>
#include "logger.h"

class Logger;

class UmbrellaParameter {
	public:
		UmbrellaParameter (char *param_vname, char *target_vname, \
				char *spring_vname, LAMMPS_NS::LAMMPS *lmp, int is_compute);
		UmbrellaParameter (const UmbrellaParameter& up);
		~UmbrellaParameter ();

		double compute_boltzmann_factor();
		double current_value;
		double last_accepted_value;
		void notify_accepted();
		void notify_accepted_debug(Logger *logger);
		void notify_rejected_debug(Logger *logger);
		char param_vname[100];
	private:
		int is_compute;
		char target_vname[100];
		char spring_vname[100];
		LAMMPS_NS::LAMMPS *lmp;

		/*
		   previous_potential is the bias potential of the most recently
		   accepted state.  current_potential is the bias potential of the most
		   recently sampled state.  (Here, "potential" === (Q - target)^2).
		*/
		//double current_potential, previous_potential;
};

#endif
