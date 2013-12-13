#ifndef PARSER_H
#define PARSER_H

#include <map>
#include <vector>
#include <mpi.h>
//#include <lammps.h>
#include <input.h>

#include "umbrella_step.h"
#include "umbrella_parameter.h"

class Parser {
	public:
		//Parser(const char *fname, LAMMPS_NS::LAMMPS *lmp);
		Parser(const char *fname, LAMMPS_NS::LAMMPS *lmp, MPI_Comm local_comm, \
				int local_rank, MPI_Comm roots_comm, int num_windows);
		~Parser();
		void parse();
		void execute_init();
		void print();
		int nsteps;
		UmbrellaStep **steps;

		int nparams;
		UmbrellaParameter **param_ptrs;

	private:
		MPI_Comm local_comm;
		int local_rank;
		MPI_Comm roots_comm;
		int num_windows;

		int process_brackets(char *line);
		std::vector<UmbrellaParameter> params;
		char fname [100];
		LAMMPS_NS::LAMMPS *lmp;
		std::vector<std::string> init_block;
		std::map<std::string, UmbrellaStep> steps_map;
};
#endif
