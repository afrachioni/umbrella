#ifndef BAROSTAT_STEP_H
#define BAROSTAT_STEP_H

#include "quantity.h"
#include "umbrella_step.h"

class BarostatStep : public UmbrellaStep {
	public:
		BarostatStep(LAMMPS_NS::LAMMPS *lmp, double probability, \
			   Quantity *temperature, char* name, Quantity *pressure);
		void execute_init();
		void execute_step();
		static double get_rate();
		static void zero_rate();
	private:
		LAMMPS_NS::bigint N;
		Quantity T, P;
		double U, V;
		static double Uold, Vold;
		static int accepted_count, count;
};
#endif
