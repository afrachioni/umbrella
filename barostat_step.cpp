#include <mpi.h>
#include <lammps.h>
#include <atom.h>
#include <library.h>

#include <stdlib.h>
#include <math.h>

#include "barostat_step.h"

double BarostatStep::Uold = 0;
double BarostatStep::Vold = 0;
int BarostatStep::accepted_count = 0;
int BarostatStep::count = 0;

BarostatStep::BarostatStep(LAMMPS_NS::LAMMPS *lmp, float probability, char* name, Global *global, double pressure) : UmbrellaStep::UmbrellaStep ( lmp, probability, name, global) {
	this->P = pressure;
	UmbrellaStep::UmbrellaStep(lmp, probability, name, global);
	is_barostat = 1;
}

void BarostatStep::execute_init() {  // XXX sets Uold, Vold once per instance
	UmbrellaStep::execute_init();
	N = static_cast<int> (lmp->atom->natoms);
	Uold = *((double *) lammps_extract_compute (lmp,(char*)"thermo_pe", 0, 0));
	Vold=*((double*)lammps_extract_variable(lmp,(char*)"lu_vol",(char*)"all"));
}

void BarostatStep::execute_step() {
	T = *((double *) lammps_extract_compute (lmp,(char*)"thermo_temp", 0, 0));
	Uold = *((double *) lammps_extract_compute (lmp,(char*)"thermo_pe", 0, 0));
	Vold = *((double *)lammps_extract_variable(lmp,(char*)"lu_vol",(char*)"all"));

	execute_block (lmp, take_step_block, global);

	U = *((double *) lammps_extract_compute (lmp,(char*)"thermo_pe", 0, 0));
	V = *((double *)lammps_extract_variable(lmp,(char*)"lu_vol",(char*)"all"));


	double exp = -(U-Uold + P*(V-Vold)/eV - N*kb*T*log(V/Vold))/(kb*T);
	double accept_rand = (double) rand() / RAND_MAX;
	int accept = log (accept_rand) < exp;
	MPI_Bcast (&accept, 1, MPI_INT, 0, global->local_comm);

	char msg[100];
	sprintf( msg, "(V - Vold)/eV: %f\tU - Uold: %f\tNkT log (V/Vold): %f",
			(V - Vold)/eV, U - Uold, N*kb*T*log(V/Vold));
	//global->debug (msg);

	//sprintf(msg, "N: %d\tk: %f\tT: %f\tlog(V/Vold): %f\n",
			//N, kb, T, log(V/Vold));
	sprintf (msg, "V: %f\tVold: %f", V, Vold);
	//global->debug(msg);
	


	if (accept) {
		if (logger != NULL) logger->comment ("accept box change");
		//global->debug("\t\t\t\t\t\taccept");
		//Uold = U;
		//Vold = V;
		execute_accept();
		++accepted_count;
	} else {
		if (logger != NULL) logger->comment ("reject box change");
		//global->debug("\t\t\t\t\t\t\treject");
		execute_reject();
	}
	++count;
}

//static
double BarostatStep::get_rate() {
	return (double)accepted_count / count;
}

//static
void BarostatStep::zero_rate() {
	accepted_count = 0;
	count = 0;
}
