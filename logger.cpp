#include "logger.h"

Logger::Logger (char *fname, int nparams, UmbrellaParameter **params) {
	strcpy (this->fname, fname);
	this->nparams = nparams;
	this->params = params;
}

Logger::~Logger() {}

void Logger::init() {
	fp = fopen (fname, "w"); // TODO catch exception
	fprintf (fp, "# Anthony's Magnificent Sampler\n");
	fprintf (fp, "#\n");
	fprintf (fp, "# Parameter names:\n");
	for (int i = 0; i < nparams; ++i)
		fprintf (fp, "#\t%-2d: %s\n", i, params[i]->param_vname);
}
	

void Logger::step_taken (int step_type, int accept) {
	fprintf (fp, "%-5d", step_type);
	for (int i = 0; i < nparams; ++i) {
		fprintf (stdout, "Logging param: %s\t%f\n", params[i]->param_vname, params[i]->current_value);
		fprintf (fp, "%-15f", params[i]->current_value);
	}
	fprintf (fp, "%-5d\n", accept);
}
