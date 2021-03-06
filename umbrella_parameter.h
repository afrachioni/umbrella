#ifndef UMBRELLA_PARAMETER_H
#define UMBRELLA_PARAMETER_H

#include <mpi.h>
#include <lammps.h>

#include "quantity.h"
#include "logger.h"

class Logger;
class Quantity;

class UmbrellaParameter {
	public:
		UmbrellaParameter (Quantity *param, Quantity *target, \
				Quantity *spring, Quantity *temp, LAMMPS_NS::LAMMPS *lmp);

		double compute_boltzmann_factor();
		void notify_accepted();
		double get_current_value();
		double get_proposed_value();
		double get_spring(), get_target();
		char *get_name();
		void notify_accepted_debug(Logger *logger);
		void notify_rejected_debug(Logger *logger);
	private:
		Quantity param, target, spring, temp;
		double current_value;
		double proposed_value;
		LAMMPS_NS::LAMMPS *lmp;
};
#endif
