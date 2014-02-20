#include <mpi.h>
#include <lammps.h>
#include <input.h> // run zero
#include <atom.h>
#include <library.h>

#include <stdlib.h>
#include <math.h>

#include "barostat_step.h"
BarostatStep::BarostatStep(LAMMPS_NS::LAMMPS *lmp, float probability, char* name, Global *global, double pressure) : UmbrellaStep::UmbrellaStep ( lmp, probability, name, global) {
	this->pressure = pressure;
	UmbrellaStep::UmbrellaStep(lmp, probability, name, global);
	N = static_cast<int> (lmp->atom->natoms);
}

void BarostatStep::execute_step() {
	global->debug("entering BarostatStep::execute_step");
	lmp->input->one("run 0");
	global->debug("Ran zero");
	T = *((double *) lammps_extract_compute (lmp,(char*)"thermo_temp", 0, 0));
	global->debug("Extracted temperature");
	Uold = *((double *) lammps_extract_compute (lmp,(char*)"thermo_pe", 0, 0));
	global->debug("Extracted U");
	P = *((double *) lammps_extract_compute (lmp,(char*)"thermo_press", 0, 0));
	global->debug("Extracted pressure");
	Vold=*((double*)lammps_extract_variable(lmp,(char*)"lu_vol",(char*)"all"));
	execute_block (lmp, take_step_block, global);
	U = *((double *) lammps_extract_compute (lmp,(char*)"thermo_pe", 0, 0));
	V = *((double *)lammps_extract_variable(lmp,(char*)"lu_vol",(char*)"all"));
	double exp = -(U-Uold + P*(V-Vold)*eV - N*kb*T*log(V/Vold))/(kb*T);
	double accept_rand = (double) rand() / RAND_MAX;
	int accept = log (accept_rand) < exp;
	MPI_Bcast (&accept, 1, MPI_INT, 0, global->local_comm);
	if (accept) {
		Uold = U; Vold = V;
		execute_accept();
	} else
		execute_reject();
}
