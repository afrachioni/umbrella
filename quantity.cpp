#include <stdlib.h> // strtod
#include <cstring>

#include <mpi.h>
#include <lammps.h>
#include <compute.h>
#include <modify.h>
#include <update.h>
#include <input.h>
#include <variable.h>

#include "quantity.h"

#include "global.h"

Quantity::Quantity(char *q, LAMMPS_NS::LAMMPS *lmp, \
		bool positive, bool integer) {
	this->lmp = lmp;
	if (q[0] == 'c' && q[1] == '_')
		compute = true;
	else if (q[0] == 'v' && q[1] == '_')
		variable = true;
	else {
		char *e;
		constant = strtod (q, &e);
		if (*e != 0) valid = false;
		if (positive && constant < 0) valid = false;
		if (integer && (long long) constant != constant) valid = false;
	}
	if (compute || variable)
		strcpy (name, q + 2);
}

/* XXX pretty sure this is default behavior
Quantity::Quantity(const Quantity &q) {
	lmp = q.lmp;
	strcpy (name, q.name);
	compute = q.compute;
	variable = q.variable;
	valid = q.valid;
	positive = q.positive;
	constant = q.constant;
}
*/

double Quantity::get_value() {
	char line[1000];
	if (!valid)
		Global::get_instance()->abort("Attempt to get value from invalid "
				"quantity!  This is an internal error which should never "
				"happen.");
	if (compute) {
		  int icompute = lmp->modify->find_compute(name);
		  if (icompute < 0) {
			  sprintf (line, "Could not find compute of name \"%s\".", name);
			  Global::get_instance()->abort(line);
		  }
		  LAMMPS_NS::Compute *compute = lmp->modify->compute[icompute];
		  if (!compute->scalar_flag || compute->local_flag || 
				  compute->peratom_flag) {
			  sprintf (line, "Compute of name \"%s\" does not compute a"
					 " global scalar.", name);
			  Global::get_instance()->abort(line);
		  }
		  if (compute->invoked_scalar != lmp->update->ntimestep)
			compute->compute_scalar();
		  return compute->scalar;
	} else if (variable) {
		int ivar = lmp->input->variable->find(name);
		if (ivar < 0) {
			sprintf (line, "Could not find variable of name \"%s\".", name);
			Global::get_instance()->abort(line);
		}
		if (!lmp->input->variable->equalstyle(ivar)) {
			sprintf (line, "Variable \"%s\" is not an equal style variable.",
					name);
			Global::get_instance()->abort(line);
		}
		return lmp->input->variable->compute_equal(ivar);
	} else
		return constant;
}

bool Quantity::is_constant() {
	return !(compute || variable) && valid;
}

bool Quantity::is_valid() {
	return valid;
}
