#include "lammps.h"
#include "quantity.h"

Quantity::Quantity(char *q, LAMMPS *lmp) {
	this->lmp = lmp;
	if (q[0] == 'c' && q[1] == '_')
		is_compute = true;
	else if (q[0] == 'v' && q[1] == '_')
		is_variable = true;
	else {
		char *e;
		constant = std::strtod (q, &e);
		if (*e != 0) valid = false;
	}
	if (is_compute || is_variable)
		strcpy (name, q + 2);
}

double Quantity::get_value() {
	if (!is_valid)
		Global::get_instance->abort("Attempt to get value from invalid "
				"quantity!  This is an internal error which should never "
				"happen.");
	if (is_compute) {
		  int icompute = lmp->modify->find_compute(name);
		  if (icompute < 0) return NULL;
		  Compute *compute = lmp->modify->compute[icompute];
		  if (!compute->scalar_flag) return NULL;
		  if (compute->invoked_scalar != lmp->update->ntimestep)
			compute->compute_scalar();
		  return compute->scalar;
	} else if (is_variable) {
		int ivar = lmp->input->variable->find(name);
		if (ivar < 0) return NULL;
		if (lmp->input->variable->equalstyle(ivar)) {
			double *dptr = (double *) malloc(sizeof(double));
			return lmp->input->variable->compute_equal(ivar);
		} else
			return NULL;
	} else
		return constant;
}

bool Quantity::is_constant() {
	return !(is_compute || is_variable) && is_valid;
}

bool Quantity::is_valid() {
	return is_valid;
}
