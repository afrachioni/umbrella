#ifndef LOGGER_H
#define LOGGER_H

#include <stdio.h>
#include <vector>

#include "umbrella_parameter.h"
#include "umbrella_step.h"


class Logger {
	public:
		Logger (char *fname, int nparams, int window_index, int local_rank, \
				UmbrellaParameter **params, int nsteps, UmbrellaStep **steps);
		~Logger ();

		void init();
		void step_taken(int step_index, int step_type, int accept);

	private:
		int nparams;
		int window_index;
		int local_rank;
		UmbrellaParameter **params;
		int nsteps;
		UmbrellaStep **steps;
		char fname[100];
		FILE *fp;
};

#endif
