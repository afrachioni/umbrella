#include <math.h> // infinity

#include "quantity.h"
#include "umbrella_parameter.h"


UmbrellaParameter::UmbrellaParameter (Quantity *param, Quantity *target, \
		Quantity *spring, Quantity *temp, LAMMPS_NS::LAMMPS *lmp) 
: param(*param), target(*target), spring(*spring), temp(*temp), lmp(lmp) {}

// ALL the physics lives here
double UmbrellaParameter::compute_boltzmann_factor() {
	proposed_value = param;
	double temperature = temp;
	if (temperature == 0) return -INFINITY; // Avoid nan
	double proposed_potential = (proposed_value-target)*(proposed_value-target);
	double current_potential = (current_value-target)*(current_value-target);
	double rval = -0.5 * spring / temperature * \
				  (proposed_potential - current_potential);
	return rval;
}

void UmbrellaParameter::notify_accepted() {
	current_value = proposed_value;
}

double UmbrellaParameter::get_current_value() {
	return current_value;
}

double UmbrellaParameter::get_proposed_value() {
	return proposed_value;
}

double UmbrellaParameter::get_spring() {
	return spring;
}

double UmbrellaParameter::get_target() {
	return target;
}

char *UmbrellaParameter::get_name() {
	return param.get_name();
}

void UmbrellaParameter::notify_accepted_debug(Logger *logger) {
	notify_accepted();

	char line[100];
	sprintf (line, "%-15s| old: %f\tnew: %f", get_name(), current_value, proposed_value);
	logger->comment (line);
}

void UmbrellaParameter::notify_rejected_debug(Logger *logger) {
	char line[100];
	sprintf (line, "%-15s| old: %f\tnew: %f", get_name(), current_value, proposed_value);
	logger->comment (line);
}
