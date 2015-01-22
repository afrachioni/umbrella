#ifndef QUANTITY_H
#define QUANTITY_H

#include "lammps.h"

class Quantity {
	public:
		Quantity(char *q, LAMMPS_NS::LAMMPS *lmp, bool positive, bool integer);
		double get_value();
		char *get_name();
		bool is_constant();
		bool is_valid();
		double operator+(Quantity &rhs);
		double operator-(Quantity &rhs);
		operator double() {return get_value();}

	private:
		LAMMPS_NS::LAMMPS *lmp;
		char name[100];
		bool compute, variable, valid;
		bool positive, integer;
		double constant;
};
#endif
