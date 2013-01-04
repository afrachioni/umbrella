#include <stdlib.h>
#include <iostream>
#include "stdio.h"
#include "optionparser.h"

#include "parser.h"

enum options_index { TARGET_MIN, TARGET_MAX, NUM_WINDOWS, START, COUNT, DURATION, TEMPERATURE, SPRING, CUTOFF, HELP };
int parser::verbose = 0;
const option::Descriptor usage [] =
{
	{ TARGET_MIN, 0, "", "target_min", parser::NumberCheck, "minimum value of order parameter"},
	{ TARGET_MAX, 0, "", "target_max", parser::NumberCheck, "maximum value of order parameter"},
	{ NUM_WINDOWS, 0, "", "num_windows", parser::IntegerCheck, "number of sampling windows" },
	{ START, 0, "", "start", parser::FileCheck, "restart file at which to begin sampling" },
	{ COUNT, 0, "", "count", parser::IntegerCheck, "number of umbrella steps" },
	{ DURATION, 0, "", "duration", parser::IntegerCheck, "duration of one umbrella step / fs" },
	{ TEMPERATURE, 0, "", "temperature", parser::NumberCheck, "temperature / K" },
	{ SPRING, 0, "", "spring", parser::NumberCheck, "spring constant / boltzmann constant"},
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
			option::printUsage(std::cerr, usage);
		parse_error ++;
		return;
	}

	int missing = 0;
	for (int i = TARGET_MIN; i < HELP; i++) {
		option::Option opt = options[i];
		if (!opt) {
			if (verbose) fprintf (stderr, "--%s is a required option!\n", usage[i].longopt);
			missing++;
		}
	}
	if (missing)
		{ parse_error ++; return; }


	if (verbose) fprintf (stderr, "Program options seem sane. Proceeding with configuration:\n");
	const char *target_min_arg = options[TARGET_MIN].last()->arg;
	target_min = atof (target_min_arg);
	if (verbose) fprintf (stderr, "\tTarget min:                     %f\n", target_min);

	const char *target_max_arg = options[TARGET_MAX].last()->arg;
	target_max = atof (target_max_arg);
	if (verbose) fprintf (stderr, "\tTarget max:                     %f\n", target_max);

	const char *num_windows_arg = options[NUM_WINDOWS].last()->arg;
	num_windows = atoi (num_windows_arg);
	if (verbose) fprintf (stderr, "\tNumber of sampling windows:     %d\n", num_windows);


	start = options[START].last()->arg;
	if (verbose) fprintf (stderr, "\tRestart filename:               %s\n", start);

	const char *count_arg = options[COUNT].last()->arg;
	count = atoi (count_arg);
	if (verbose) fprintf (stderr, "\tNumber of umbrella steps:       %d\n", count);

	const char *duration_arg = options[DURATION].last()->arg;
	duration = atoi (duration_arg);
	if (verbose) fprintf (stderr, "\tUmbrella step duration:         %d\n", duration);

	const char *temperature_arg = options[TEMPERATURE].last()->arg;
	temperature = atof (temperature_arg);
	if (verbose) fprintf (stderr, "\tTemperature:                    %f\n", temperature);

	const char *spring_arg = options[SPRING].last()->arg;
	spring = atof (spring_arg);
	if (verbose) fprintf (stderr, "\tSpring constant:                %f\n", spring);

	const char *cutoff_arg = options[CUTOFF].last()->arg;
	cutoff = atof (cutoff_arg);
	if (verbose) fprintf (stderr, "\tQ6 neighbor cutoff:             %f\n", cutoff);
}
option::ArgStatus parser::MyCheck (const option::Option& option, bool msg)
{
	if (verbose) fprintf (stderr, "option->arg: %s\n", option.arg);
	return option::ARG_ILLEGAL;

}
option::ArgStatus parser::FileCheck (const option::Option& option, bool msg)
{
	const char *arg = option.arg;
	if (arg == 0 || arg[0] == '\0') {
		if (msg && verbose) fprintf (stderr, "Option '%s' requires an argument!\n", option.name);
		return option::ARG_ILLEGAL;
	}
	FILE *f = fopen (arg, "r");
	if (f == NULL) {
		if (msg && verbose)
			fprintf (stderr, "Unable to read file: %s\n", arg);
		return option::ARG_ILLEGAL;
	} else {
		return option::ARG_OK;
	}
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
