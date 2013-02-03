#ifndef PARSER_H
#define PARSER_H

#include "optionparser.h"
#include "stdio.h"
#define WHAM_BOLTZMANN 0.0083144621 // kJ/mol-K 

class parser
{
	public:
		parser (int narg, char **arg);
		static int verbose;
		int parse_error;
		FILE *init;
		int count;
		int duration;
		double temperature;
		double spring;
		int l;
		double cutoff;

		static option::ArgStatus MyCheck (const option::Option& option, bool msg);
		static option::ArgStatus NumberCheck (const option::Option& option, bool msg);
		static option::ArgStatus IntegerCheck (const option::Option& option, bool msg);
		static option::ArgStatus NonEmpty (const option::Option& option, bool msg);
};
#endif
