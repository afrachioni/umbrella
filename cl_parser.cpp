#include <stdlib.h>
#include <iostream>
#include <string.h>
#include "stdio.h"
#include "optionparser.h"

#include "cl_parser.h"

enum options_index { COUNT, WINDOWS, LOG_LAMMPS, SCRIPT, HELP };
int CLParser::verbose = 0;
//char CLParser::err_str[1000];
const option::Descriptor usage [] =
{
	{ COUNT, 0, "c", "count", CLParser::IntegerCheck, "\t--count\t     number of umbrella steps" },
	{ WINDOWS, 0, "w", "windows", CLParser::IntegerCheck, "\t--windows\t     number of parallel windows"},
	{ LOG_LAMMPS, 0, "l", "lammps_logs", option::Arg::None, "\t--lammps_log\t     log LAMMPS IO to log_*.lammps"},
	{ SCRIPT, 0, "f", "script", CLParser::NonEmpty, "\t--script\t     name of input script"},
	{ HELP, 0, "", "help", option::Arg::None, "\t--help\t     print usage" },
	{ 0, 0, 0, 0, 0, 0 }
};
CLParser::CLParser (int narg, char **arg)
{
	parse_error = 0;
	//strcpy (err_str, "");
	//char line[100];
	narg-=(narg>0); arg+=(narg>0); // skip program name arg[0] if present
	option::Stats  stats(usage, narg, arg);
	option::Option options[stats.options_max], buffer[stats.buffer_max];
	option::Parser parse(usage, narg, arg, options, buffer);

	if (parse.error()) // Who calls this?
		{++parse_error; return;}

	if (options[HELP] || narg == 0) {
		if (verbose)
			option::printUsage(std::cout, usage);
		++parse_error;
		return;
	}

	if (options[LOG_LAMMPS])
		log_lammps = 1;

	int missing = 0;
	for (int i = 0; i < HELP; i++) {
		option::Option opt = options[i];
		if (!opt && i != LOG_LAMMPS) { //XXX a cleaner way probably exists
			//sprintf ("--%s is a required option!\n", usage[i].longopt);
			//strcat (err_str, line);
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

	const char *windows_arg = options[WINDOWS].last()->arg;
	windows = atoi (windows_arg);
	if (verbose) fprintf (stdout, "\tNumber of parallel windows:            %d\n", windows);

	strcpy (script, options[SCRIPT].last()->arg);
	FILE *init = fopen (script, "r");
	if (init == NULL) {
		if (verbose) fprintf (stderr, "Unable to read input script: %s\n",  script);
		++parse_error;
		return;
	}
	if (verbose) fprintf (stdout, "\tInput script:                          %s\n", script);
	fclose (init);

}
// Necessary?
option::ArgStatus CLParser::MyCheck (const option::Option& option, bool msg)
{
	if (verbose) fprintf (stderr, "option->arg: %s\n", option.arg);
	return option::ARG_ILLEGAL;

}

// Unused at the moment
option::ArgStatus CLParser::NumberCheck (const option::Option& option, bool msg)
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
option::ArgStatus CLParser::IntegerCheck (const option::Option& option, bool msg)
{
	const char *arg = option.arg;
	if (arg == 0 || arg[0] == '\0') {
		char line[100];
		//sprintf (line, "Option \"%s\" requires an argument!\n");
		//strcat (err_msg, line);
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
option::ArgStatus CLParser::NonEmpty(const option::Option& option, bool msg)
{
	if (option.arg != 0 && option.arg[0] != 0)
		return option::ARG_OK;

	if (msg) fprintf (stderr, "Option '%s' requires a non-empty argument\n", option.name);
	return option::ARG_ILLEGAL;
}
