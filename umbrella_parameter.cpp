#include <mpi.h>
#include <lammps.h>
#include <library.h>
#include <math.h> // infinity
#include <string.h>
#include "quantity.h"
#include "umbrella_parameter.h"


UmbrellaParameter::UmbrellaParameter (Quantity *param, Quantity *target, \
		Quantity *spring, LAMMPS_NS::LAMMPS *lmp) {
	this->param_Q = new Quantity(*param);
	this->target_Q = new Quantity(*target);
	this->spring_Q = new Quantity(*spring);
	this->lmp = lmp;
}

/*
// XXX this might be default behavior anyway
UmbrellaParameter::UmbrellaParameter (const UmbrellaParameter& up) {
	//this->param_Q = up.param_Q;
	//this->target_Q = up.target_Q;
	//this->spring_Q = up.spring_Q;
	this->lmp = up.lmp;
}
*/

// ALL the physics lives here
double UmbrellaParameter::compute_boltzmann_factor() {
	current_value = param_Q->get_value();

	double temperature = *((double *) lammps_extract_compute(lmp,(char*)"thermo_temp", 0, 0));
	//temperature = 1; //XXX
	//temperature = 1/8.617e-5; //XXX
	if (temperature == 0) return -INFINITY; // Avoid nan

	double target = target_Q->get_value();
	double spring = spring_Q->get_value();
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

double UmbrellaParameter::get_spring() {
	return spring_Q->get_value();
}

double UmbrellaParameter::get_target() {
	return target_Q->get_value();
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

UmbrellaParameter::~UmbrellaParameter() {};
