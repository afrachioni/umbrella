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
		int N; // XXX shouldn't this be an lmptype or something?
		Quantity *pressure;
		double T, U, P, V;
		static double Uold, Vold;
		static int accepted_count, count;
};
#endif
