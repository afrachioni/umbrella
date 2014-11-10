#ifndef PARSER_H
#define PARSER_H

#include <map>
#include <vector>
#include <mpi.h>
//#include <lammps.h>
#include <input.h>

#include "global.h"
#include "quantity.h"
#include "umbrella_step.h"
#include "umbrella_parameter.h"
#include "periodic_task.h"
#include "histogram.h"

class Parser {
	public:
		Parser(const char *fname, LAMMPS_NS::LAMMPS *lmp, Global *global);
		~Parser();
		int parse();
		void execute_init();
		void print();
		int nsteps;
		UmbrellaStep **steps;
		double get_temp();
		char *error_message();

		int nparams;
		UmbrellaParameter **param_ptrs;
		std::vector<PeriodicTask*> tasks; // XXX these can be bare objects, I think
		std::vector<Histogram*> histograms;

		int bias_every;

		std::map<std::string, UmbrellaStep*> steps_map;
	private:
		char msg[100];
		void process_brackets(char *line);
		std::vector<UmbrellaParameter *> params;
		char fname [100];
		Quantity *temp;
		LAMMPS_NS::LAMMPS *lmp;
		Global *global;
		std::vector<std::string> init_block;
};
#endif
