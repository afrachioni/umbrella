#ifndef UMBRELLA_STEP_H
#define UMBRELLA_STEP_H

#include <mpi.h>
#include <lammps.h>
#include <string>
#include<vector>

#include "quantity.h"
#include "global.h"
#include "logger.h"

class Logger;

class UmbrellaStep {
	public:
		int is_barostat;

		UmbrellaStep(LAMMPS_NS::LAMMPS *lmp, Quantity *probability, char *name);
		UmbrellaStep();
		~UmbrellaStep();
		std::vector<std::string>* get_step_init_block();
		std::vector<std::string>* get_take_step_block();
		std::vector<std::string>* get_if_accept_block();
		std::vector<std::string>* get_if_reject_block();

		virtual void execute_init();
		virtual void execute_step();
		void execute_accept();
		void execute_reject();

		Quantity *probability;
		float rand_min;
		float rand_max;

		char name[100];
		static int get_atoms_called;
		static int get_types_called;
		static int force_accept;
		static void execute_block(LAMMPS_NS::LAMMPS *lmp, std::vector<std::string> block);
//ex-private (needed for derived BarostatStep
		LAMMPS_NS::LAMMPS *lmp;

		static double *positions_buffer;
		static int *types_buffer;
		static double xlo, xhi, ylo, yhi, zlo, zhi;
		static char line[];

		std::vector<std::string> step_init_block;
		std::vector<std::string> take_step_block;
		std::vector<std::string> if_accept_block;
		std::vector<std::string> if_reject_block;

		Logger *logger;
		void set_logger_debug (Logger *logger);
};
#endif
