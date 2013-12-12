#include "logger.h"

Logger::Logger (char *fname, int nparams, UmbrellaParameter **params,\
		int nsteps, UmbrellaStep **steps) {
	strcpy (this->fname, fname);
	this->nparams = nparams;
	this->params = params;
	this->nsteps = nsteps;
	this->steps = steps;
}

Logger::~Logger() {}

void Logger::init() {
	fp = fopen (fname, "w"); // TODO catch exception
	fprintf (fp, "# Anthony's Magnificent Sampler\n");
	fprintf (fp, "#\n");
	fprintf (fp, "# Step types:\n");
	for (int i = 0; i < nsteps; ++i)
		fprintf (fp, "#\t%-2d: %s\n", i, steps[i]->name);
	fprintf (fp, "#\n");


	fprintf (fp, "#Step ");
	for (int i = 0; i < nparams; ++i)
		fprintf (fp, "%-15s", params[i]->param_vname);
	fprintf (fp, "Type  Accept\n");
}
	

void Logger::step_taken (int step_index, int step_type, int accept) {
	fprintf (fp, "%-6d", step_index);
	for (int i = 0; i < nparams; ++i) {
		fprintf (fp, "%-15f", params[i]->current_value);
	}
	fprintf (fp, "%-5d", step_type);
	fprintf (fp, "%-7d\n", accept);
}
