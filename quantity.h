#ifndef QUANTITY_H
#define QUANTITY_H

#include "lammps.h"

class Quantity {
	public:
		Quantity(char *q, LAMMPS_NS::LAMMPS *lmp);
		double get_value();
		bool is_constant();
		bool is_valid();

	private:
		LAMMPS_NS::LAMMPS *lmp;
		char name[100];
		bool compute, variable, valid;
		double constant;
};
#endif
