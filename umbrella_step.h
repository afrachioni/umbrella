#ifndef UMBRELLA_STEP_H
#define UMBRELLA_STEP_H

#include <mpi.h>
#include <lammps.h>

#include "global.h"


class UmbrellaStep {
	public:
		UmbrellaStep(LAMMPS_NS::LAMMPS *lmp, float probability, char *name, Global *global);
		UmbrellaStep();
		~UmbrellaStep();
		std::vector<std::string>* get_step_init_block();
		std::vector<std::string>* get_take_step_block();
		std::vector<std::string>* get_if_accept_block();
		std::vector<std::string>* get_if_reject_block();

		void execute_init();
		void execute_step();
		void execute_accept();
		void execute_reject();

		float probability;
		float rand_min;
		float rand_max;

		char name[100];
		static int get_atoms_called;
		static int force_accept;
		static void execute_block(LAMMPS_NS::LAMMPS *lmp, std::vector<std::string> block, Global *global);
	private:
		LAMMPS_NS::LAMMPS *lmp;
		static double *positions_buffer;
		Global *global;
		std::vector<std::string> step_init_block;
		std::vector<std::string> take_step_block;
		std::vector<std::string> if_accept_block;
		std::vector<std::string> if_reject_block;
};
#endif
