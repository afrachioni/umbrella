#include <mpi.h>
#include <lammps.h>
#include <atom.h>
#include <library.h>

#include <stdlib.h>
#include <math.h>

#include "quantity.h"
#include "barostat_step.h"

double BarostatStep::Uold = 0;
double BarostatStep::Vold = 0;
int BarostatStep::accepted_count = 0;
int BarostatStep::count = 0;

BarostatStep::BarostatStep(LAMMPS_NS::LAMMPS *lmp, Quantity *probability, \
		char* name, Quantity *pressure) :\
		UmbrellaStep::UmbrellaStep (lmp, probability, name) {
	this->pressure = new Quantity(*pressure);
	is_barostat = 1;
}

void BarostatStep::execute_init() {  // XXX sets Uold, Vold once per instance
	UmbrellaStep::execute_init();
	N = static_cast<int> (lmp->atom->natoms);
	Uold = *((double *) lammps_extract_compute (lmp,(char*)"thermo_pe", 0, 0));
	Vold=*((double*)lammps_extract_variable(lmp,(char*)"lu_vol",(char*)"all"));
}

void BarostatStep::execute_step() {
	// TODO might be able to store pointer, dereference when I need it, rather
	// than calling extract each step (not sure whether this matters)
	T = *((double *) lammps_extract_compute (lmp, (char*)"thermo_temp", 0, 0));
	Uold = *((double *) lammps_extract_compute (lmp,(char*)"thermo_pe", 0, 0));
	Vold=*((double*)lammps_extract_variable(lmp,(char*)"lu_vol",(char*)"all"));

	execute_block (lmp, take_step_block);

	U = *((double *) lammps_extract_compute (lmp,(char*)"thermo_pe", 0, 0));
	V = *((double *)lammps_extract_variable(lmp,(char*)"lu_vol",(char*)"all"));


	// XXX This is the long range correction to the pressure, see Verlet
	//double density = *((double *) lammps_extract_variable(lmp,(char*)"density",(char*)"all"));
	//double density_factor = density*density * 1.0694;
	//double density_factor = density * 0.929915;
	double density_factor = 0;

	P = pressure->get_value();
	double exp = -(U-Uold + (P + density_factor)*(V-Vold)/eV - N*kb*T*log(V/Vold))/(kb*T);
	double accept_rand = (double) rand() / RAND_MAX;
	int accept = log (accept_rand) < exp;
	MPI_Bcast (&accept, 1, MPI_INT, 0, Global::get_instance()->local_comm);

	if (accept) {
		if (logger != NULL) logger->comment ((char*)"accept box change");
		execute_accept();
		++accepted_count;
	} else {
		if (logger != NULL) logger->comment ((char*)"reject box change");
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
