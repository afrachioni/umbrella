#include <mpi.h>
#include <lammps.h>
#include <library.h>
#include <math.h> // infinity
#include <string.h>
#include "umbrella_parameter.h"


UmbrellaParameter::UmbrellaParameter (char *param_vname, char *target_vname, \
		char *spring_vname, LAMMPS_NS::LAMMPS *lmp, int is_compute) {
	strcpy (this->param_vname, param_vname);
	strcpy (this->target_vname, target_vname);
	strcpy (this->spring_vname, spring_vname);
	this->lmp = lmp;
	this->is_compute = is_compute;
}

UmbrellaParameter::UmbrellaParameter (const UmbrellaParameter& up) {
	strcpy (this->param_vname, up.param_vname);
	strcpy (this->target_vname, up.target_vname);
	strcpy (this->spring_vname, up.spring_vname);
	this->lmp = up.lmp;
	this->is_compute = up.is_compute;
}

// ALL the physics lives here
double UmbrellaParameter::compute_boltzmann_factor() {
	if (is_compute)
		current_value = *((double *) lammps_extract_compute(lmp,param_vname, 0, 0));
	else
		current_value = *((double *) lammps_extract_variable(lmp, \
					param_vname, (char *) "all")); //TODO pass group in

	double temperature = *((double *) lammps_extract_compute(lmp,(char*)"thermo_temp", 0, 0));
	//temperature = 1; //XXX
	//temperature = 1/8.617e-5; //XXX
	if (temperature == 0) return -INFINITY; // Avoid nan
	
	/*
	//double target = *((double *) lammps_extract_variable(lmp, \
				target_vname, (char *) "all")); //TODO pass group in
	//double spring = *((double *) lammps_extract_variable(lmp, \
				spring_vname, (char *) "all")); //TODO pass group in
	*/
	double target = extract_target();
	double spring = extract_spring();
	double current_potential = (current_value-target)*(current_value-target);
	double previous_potential = (last_accepted_value-target)*(last_accepted_value-target);
	double rval = -0.5 * spring / temperature * \
				  (current_potential - previous_potential);
	return rval;
}

void UmbrellaParameter::notify_accepted() {
	last_accepted_value = current_value;
}

double UmbrellaParameter::get_current_value() {
	return current_value;
}

double UmbrellaParameter::get_last_accepted_value() {
	return last_accepted_value;
}

void UmbrellaParameter::notify_accepted_debug(Logger *logger) {
	notify_accepted();

	char line[100];
	sprintf (line, "%-15s| old: %f\tnew: %f", param_vname, last_accepted_value, current_value);
	logger->comment (line);
}

void UmbrellaParameter::notify_rejected_debug(Logger *logger) {
	char line[100];
	sprintf (line, "%-15s| old: %f\tnew: %f", param_vname, last_accepted_value, current_value);
	logger->comment (line);
}

double UmbrellaParameter::extract_target() {
	return *((double *) lammps_extract_variable(lmp, \
				target_vname, (char *) "all")); //TODO pass group in
}

double UmbrellaParameter::extract_spring() {
	return *((double *) lammps_extract_variable(lmp, \
				spring_vname, (char *) "all")); //TODO pass group in
}

UmbrellaParameter::~UmbrellaParameter() {};
