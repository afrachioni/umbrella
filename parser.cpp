#include <stdlib.h>
#include <iostream>
#include "stdio.h"
#include "optionparser.h"

#include "parser.h"

enum options_index { INIT, COUNT, DURATION, TEMPERATURE, L, CUTOFF, HELP };
int parser::verbose = 0;
const option::Descriptor usage [] =
{
	{ INIT, 0, "", "init", parser::NonEmpty, "file containing paths to initial restart files and window centers" },
	{ COUNT, 0, "", "count", parser::IntegerCheck, "number of umbrella steps" },
	{ DURATION, 0, "", "duration", parser::IntegerCheck, "duration of one umbrella step / fs" },
	{ TEMPERATURE, 0, "", "temperature", parser::NumberCheck, "temperature / K" },
	{ L, 0, "", "l", parser::IntegerCheck, "l for order parameter Q" },
	{ CUTOFF, 0, "", "cutoff", parser::NumberCheck, "neighbor cutoff for Q6" },
	{ HELP, 0, "", "help", option::Arg::None, "help usage" },
	{ 0, 0, 0, 0, 0, 0 }
};
parser::parser (int narg, char **arg)
{
	parse_error = 0;
	narg-=(narg>0); arg+=(narg>0); // skip program name arg[0] if present
	option::Stats  stats(usage, narg, arg);
	option::Option options[stats.options_max], buffer[stats.buffer_max];
	option::Parser parse(usage, narg, arg, options, buffer);

	if (parse.error())
		{parse_error ++; return;}

	if (options[HELP] || narg == 0) {
		if (verbose)
			option::printUsage(std::cout, usage);
		parse_error ++;
		return;
	}

	int missing = 0;
	for (int i = 0; i < HELP; i++) {
		option::Option opt = options[i];
		if (!opt) {
			if (verbose) fprintf (stderr, "--%s is a required option!\n", usage[i].longopt);
			missing++;
		}
	}
	if (missing)
		{ ++parse_error; return; }


	if (verbose) fprintf (stdout, "Program options seem sane. Proceeding with configuration:\n");

	const char *count_arg = options[COUNT].last()->arg;
	count = atoi (count_arg);
	if (verbose) fprintf (stdout, "\tNumber of umbrella steps:              %d\n", count);

	const char *duration_arg = options[DURATION].last()->arg;
	duration = atoi (duration_arg);
	if (verbose) fprintf (stdout, "\tUmbrella step duration:                %d\n", duration);

	const char *temperature_arg = options[TEMPERATURE].last()->arg;
	temperature = atof (temperature_arg);
	if (verbose) fprintf (stdout, "\tTemperature:                           %f\n", temperature);

	const char *l_arg = options[L].last()->arg;
	l = atoi (l_arg);
	if (verbose) fprintf (stdout, "\tl:                                     %d\n", l);

	const char *cutoff_arg = options[CUTOFF].last()->arg;
	cutoff = atof (cutoff_arg);
	if (verbose) fprintf (stdout, "\tQ6 neighbor cutoff:                    %f\n", cutoff);

	const char *init_arg = options[INIT].last()->arg;
	init = fopen (init_arg, "r");
	if (init == NULL) {
		//TODO this line prints multiple times for some reason
		fprintf (stderr, "Unable to read file: %s\n", init_arg); ++parse_error; return;
	}
	if (verbose) fprintf (stdout, "\tInitial window configuration file:     %s\n", init_arg);
}
option::ArgStatus parser::MyCheck (const option::Option& option, bool msg)
{
	if (verbose) fprintf (stderr, "option->arg: %s\n", option.arg);
	return option::ARG_ILLEGAL;

}
option::ArgStatus parser::NumberCheck (const option::Option& option, bool msg)
{
	const char *arg = option.arg;
	if (arg == 0 || arg[0] == '\0') {
		if (msg) fprintf (stderr, "Option '%s' requires an argument!\n", option.name);
		return option::ARG_ILLEGAL;
	}
	int periods = 0;
	int err = 0;
	char c = 'e';
	for (int i = 0; arg[i] != '\0'; i++) {
		if (arg[i] == '.')
			periods++;
		else if (periods > 1 || arg[i] < '0' || arg[i] > '9')
			err++;
	}
	if (err) {
		if (msg && verbose)
			fprintf (stderr, "Argument to '%s' doesn't look like a number: %s\n", option.name, arg);
		return option::ARG_ILLEGAL;
	} else
		return option::ARG_OK;
}
option::ArgStatus parser::IntegerCheck (const option::Option& option, bool msg)
{
	const char *arg = option.arg;
	if (arg == 0 || arg[0] == '\0') {
		if (msg && verbose) fprintf (stderr, "Option '%s' requires an argument!\n", option.name);
		return option::ARG_ILLEGAL;
	}
	int periods = 0;
	int err = 0;
	char c = 'e';
	for (int i = 0; arg[i] != '\0'; i++) {
		if (arg[i] == '.')
			periods++;
		else if (periods > 0 || arg[i] < '0' || arg[i] > '9')
			err++;
	}
	if (err) {
		if (msg && verbose)
			fprintf (stderr, "Argument to '%s' doesn't look like an integer: %s\n", option.name, arg);
		return option::ARG_ILLEGAL;
	} else
		return option::ARG_OK;
}
option::ArgStatus parser::NonEmpty(const option::Option& option, bool msg)
{
	if (option.arg != 0 && option.arg[0] != 0)
		return option::ARG_OK;

	if (msg) fprintf (stderr, "Option '%s' requires a non-empty argument\n", option.name);
	return option::ARG_ILLEGAL;
}
