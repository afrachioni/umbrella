#ifndef BAROSTAT_STEP_H
#define BAROSTAT_STEP_H

#include "umbrella_step.h"

class BarostatStep : public UmbrellaStep {
	public:
		BarostatStep(LAMMPS_NS::LAMMPS *lmp, float probability, char* name, Global *global, double pressure);
		void execute_step();
	private:
		double T, U, P, V, N, Uold, Vold, pressure;
		static const double kb = 8.61733238e-5; // eV / K
		static const double eV = 1.6021766e6; // bar A^3

		// From UmbrellaStep
		//LAMMPS_NS::LAMMPS *lmp;
		//std::vector<std::string> take_step_block;
		//Global *global;
};
#endif
