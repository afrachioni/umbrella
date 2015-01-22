#include <mpi.h>
#include <lammps.h>
#include <library.h>
#include <math.h> // infinity
#include <string.h>
#include "quantity.h"
#include "umbrella_parameter.h"


UmbrellaParameter::UmbrellaParameter (Quantity *param, Quantity *target, \
		Quantity *spring, Quantity *temp, LAMMPS_NS::LAMMPS *lmp) 
: param(*param), target(*target), spring(*spring), temp(*temp), lmp(lmp) {}

// ALL the physics lives here
double UmbrellaParameter::compute_boltzmann_factor() {
	current_value = param.get_value();
	double temperature = temp.get_value();
	if (temperature == 0) return -INFINITY; // Avoid nan
	double current_potential = (current_value-target)*(current_value-target);
	double previous_potential = (last_accepted_value - target) * \
								(last_accepted_value - target);
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
	return spring.get_value();
}

double UmbrellaParameter::get_target() {
	return target.get_value();
}

char *UmbrellaParameter::get_name() {
	return param.get_name();
}

void UmbrellaParameter::notify_accepted_debug(Logger *logger) {
	notify_accepted();

	char line[100];
	sprintf (line, "%-15s| old: %f\tnew: %f", get_name(), last_accepted_value, current_value);
	logger->comment (line);
}

void UmbrellaParameter::notify_rejected_debug(Logger *logger) {
	char line[100];
	sprintf (line, "%-15s| old: %f\tnew: %f", get_name(), last_accepted_value, current_value);
	logger->comment (line);
}
