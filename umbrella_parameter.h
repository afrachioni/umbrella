#ifndef UMBRELLA_PARAMETER_H
#define UMBRELLA_PARAMETER_H

#include <mpi.h>
#include <lammps.h>
#include "quantity.h"
#include "logger.h"

class Logger;

class UmbrellaParameter {
	public:
		UmbrellaParameter (Quantity *param, Quantity *target, \
				Quantity *spring, LAMMPS_NS::LAMMPS *lmp);

		double compute_boltzmann_factor();
		void notify_accepted();
		double get_current_value();
		double get_last_accepted_value();
		double get_spring(), get_target();
		char *get_name();
		void notify_accepted_debug(Logger *logger);
		void notify_rejected_debug(Logger *logger);
	private:
		Quantity *param_Q, *target_Q, *spring_Q;
		double current_value;
		double last_accepted_value;
		LAMMPS_NS::LAMMPS *lmp;
};

#endif
