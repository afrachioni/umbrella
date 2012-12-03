#ifndef PARSER_H
#define PARSER_H

#include "optionparser.h"
#include "stdio.h"

class parser
{
	public:
		parser (int narg, char **arg);
		static int verbose;
		int parse_error;
		double target_min;
		double target_max;
		int num_windows;
		const char *start;
		int count;
		int duration;
		double temperature;
		double spring;
		double cutoff;

		static option::ArgStatus MyCheck (const option::Option& option, bool msg);
		static option::ArgStatus FileCheck (const option::Option& option, bool msg);
		static option::ArgStatus NumberCheck (const option::Option& option, bool msg);
		static option::ArgStatus IntegerCheck (const option::Option& option, bool msg);
};
//int parser::verbose;
#endif
