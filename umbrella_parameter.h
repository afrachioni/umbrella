#ifndef UMBRELLA_PARAMETER_H
#define UMBRELLA_PARAMETER_H

#include <mpi.h>
#include <lammps.h>

class UmbrellaParameter {
	public:
		UmbrellaParameter (char *param_vname, char *target_vname, \
				char *spring_vname, LAMMPS_NS::LAMMPS *lmp, int is_compute);
		UmbrellaParameter (const UmbrellaParameter& up);
		~UmbrellaParameter ();

		double compute_boltzmann_factor();
		double current_value;
		char param_vname[100];
	private:
		int is_compute;
		char target_vname[100];
		char spring_vname[100];
		LAMMPS_NS::LAMMPS *lmp;

		double previous_potential;
};

#endif
