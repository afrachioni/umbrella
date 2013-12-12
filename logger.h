#ifndef LOGGER_H
#define LOGGER_H

#include <stdio.h>
#include <vector>

#include "umbrella_parameter.h"


class Logger {
	public:
		Logger (char *fname, int nparams, UmbrellaParameter **params);
		~Logger ();

		void init();
		void step_taken(int step_type, int accept);

	private:
		int nparams;
		UmbrellaParameter **params;
		char fname[100];
		FILE *fp;
};

#endif
