#include <inttypes.h>
#include <sys/time.h>
#include <string.h>
#include "logger.h"

Logger::Logger (char *fname, int nparams, int window_index, int local_rank, \
		UmbrellaParameter **params, int nsteps, UmbrellaStep **steps) {
	strcpy (this->fname, fname);
	this->nparams = nparams;
	this->window_index = window_index;
	this->local_rank = local_rank;
	this->params = params;
	this->nsteps = nsteps;
	this->steps = steps;
	this->init_time = get_time();
}

Logger::~Logger() {}

void Logger::init() {
	if (local_rank == 0) {
		fp = fopen (fname, "w"); // TODO catch exception
		fprintf (fp, "# Anthony's Magnificent Sampler\n");
		fprintf (fp, "# Window index: %d\n", window_index);
		fprintf (fp, "#\n");
		fprintf (fp, "# Step label   Name           Probabilitiy\n");
		for (int i = 0; i < nsteps; ++i)
			fprintf (fp, "#     %-2d       %-15s%f\n", i, steps[i]->name, steps[i]->probability);
		fprintf (fp, "#\n");
		fprintf (fp, "#Step    ");
		for (int i = 0; i < nparams; ++i) {
			if (i < 2) {
				char last_accepted_name[100];
				sprintf (last_accepted_name, "Last %s", params[i]->param_vname);
				fprintf (fp, "%-15s", last_accepted_name);
			}
			fprintf (fp, "%-15s", params[i]->param_vname);
		}
		fprintf (fp, "Type  Accept  Walltime\n");
	}
}

void Logger::step_taken (int step_index, int step_type, int accept) {
	if (local_rank == 0) {
		fprintf (fp, "%-9d", step_index);
		for (int i = 0; i < nparams; ++i) {
			if (i < 2)  // XXX Debug!
				fprintf(fp, "%-15f", params[i]->last_accepted_value);
			fprintf (fp, "%-15f", params[i]->current_value);
		}
		//fprintf (fp, "%-6d", step_type);
		fprintf (fp, "%-8d", accept);
		fprintf (fp, "%-10ld\n", get_time() - init_time);
		//TODO flush periodic in time, not samples
		// DEBUG if (step_index % 100 == 0)
			fflush (fp);
	}
}

void Logger::comment (char *line) {
	if (local_rank == 0)
		fprintf (fp, "# %s\n", line);
}

int64_t Logger::get_time() {
	struct timeval tv;
	gettimeofday(&tv, NULL);
	uint64_t ret = tv.tv_usec;
	ret /= 1000;
	ret += (tv.tv_sec * 1000);
	return ret;
}
