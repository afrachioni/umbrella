#ifndef QUANTITY_H
#define QUANTITY_H

#include "lammps.h"

class Quantity {
	public:
		Quantity::Quantity(char *q, LAMMPS *lmp);
		double Quantity::get_value();
		bool Quantity::is_constant();
		bool Quantity::is_valid();

	private:
		LAMMPS *lmp;
		char name[100];
		bool is_compute, is_variable, is_valid;
		double constant;
};
#endif
