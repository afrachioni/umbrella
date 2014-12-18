#ifndef BAROSTAT_STEP_H
#define BAROSTAT_STEP_H

#include "quantity.h"
#include "umbrella_step.h"

class BarostatStep : public UmbrellaStep {
	public:
		BarostatStep(LAMMPS_NS::LAMMPS *lmp, Quantity *probability, char* name, Quantity *pressure);
		void execute_init();
		void execute_step();
		static double get_rate();
		static void zero_rate();
	private:
		int N;
		Quantity *pressure;
		double T, U, P, V;
		static double Uold, Vold;
		static int accepted_count, count;

		// Metal units
		//static const double kb = 8.61733238e-5; // eV / K
		//static const double eV = 1.6021766e6; // bar A^3

		// Reduced LJ units
		static const double kb = 1;
		static const double eV = 1;
};
#endif
