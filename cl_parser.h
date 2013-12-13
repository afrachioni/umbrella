#ifndef CL_PARSER_H
#define CL_PARSER_H

#include "optionparser.h"
#include "stdio.h"

class CLParser
{
	public:
		CLParser (int narg, char **arg);
		static int verbose;
		int parse_error;
		int count;
		int windows;
		char script[100];

		static option::ArgStatus MyCheck (const option::Option& option, bool msg);
		static option::ArgStatus NumberCheck (const option::Option& option, bool msg);
		static option::ArgStatus IntegerCheck (const option::Option& option, bool msg);
		static option::ArgStatus NonEmpty (const option::Option& option, bool msg);
};
#endif
