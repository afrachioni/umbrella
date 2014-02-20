#ifndef PARSER_H
#define PARSER_H

#include <map>
#include <vector>
#include <mpi.h>
//#include <lammps.h>
#include <input.h>

#include "global.h"
#include "umbrella_step.h"
#include "umbrella_parameter.h"
#include "periodic_task.h"

class Parser {
	public:
		Parser(const char *fname, LAMMPS_NS::LAMMPS *lmp, Global *global);
		~Parser();
		void parse();
		void execute_init();
		void print();
		int nsteps;
		UmbrellaStep **steps;

		int nparams;
		UmbrellaParameter **param_ptrs;
		std::vector<PeriodicTask*> tasks;

	private:
		int process_brackets(char *line);
		std::vector<UmbrellaParameter> params;
		char fname [100];
		LAMMPS_NS::LAMMPS *lmp;
		Global *global;
		std::vector<std::string> init_block;
		std::map<std::string, UmbrellaStep*> steps_map;
};
#endif
