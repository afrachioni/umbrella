#ifndef CL_PARSER_H
#define CL_PARSER_H

#include "optionparser.h"
#include "stdio.h"

//#define WHAM_BOLTZMANN 0.001982923700 // Boltzmann's constant in kcal/mol K
#define WHAM_BOLTZMANN 0.0083144621 // kJ/mol-K 

class CLParser
{
	public:
		CLParser (int narg, char **arg);
		static int verbose;
		int parse_error;
		FILE *init;
		int count;
		int duration;
		double temperature;
		int l;
		double cutoff;

		static option::ArgStatus MyCheck (const option::Option& option, bool msg);
		static option::ArgStatus NumberCheck (const option::Option& option, bool msg);
		static option::ArgStatus IntegerCheck (const option::Option& option, bool msg);
		static option::ArgStatus NonEmpty (const option::Option& option, bool msg);
};
#endif