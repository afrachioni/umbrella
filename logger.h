#ifndef LOGGER_H
#define LOGGER_H

#include <inttypes.h>

#include "umbrella_parameter.h"
#include "umbrella_step.h"

class UmbrellaParameter;
class UmbrellaStep;

class Logger {
	public:
		Logger (char *fname, int nparams, int window_index, int local_rank, \
				UmbrellaParameter **params, int nsteps, UmbrellaStep **steps);
		~Logger ();

		void init();
		void step_taken(int step_index, int step_type, int accept);
		void comment (char *line);

		static int64_t get_time(); // could be static
	private:
		int nparams;
		int window_index;
		int local_rank;
		UmbrellaParameter **params;
		int nsteps;
		UmbrellaStep **steps;
		char fname[100];
		FILE *fp;
		int64_t init_time;
};

#endif
