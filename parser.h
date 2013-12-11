#include <map>
#include <vector>
#include <mpi.h>
//#include <lammps.h>
#include <input.h>

#include "umbrella_step.h"
#include "umbrella_parameter.h"

class Parser {
	public:
		//Parser(const char *fname, LAMMPS_NS::LAMMPS *lmp);
		Parser(const char *fname, void *lmp);
		~Parser();
		void parse();
		void execute_init();
		void print();
	private:
		char fname [100];
		LAMMPS_NS::LAMMPS *lmp;
		std::vector<std::string> init_block;
		std::map<std::string, UmbrellaStep> steps_map;
		std::vector<UmbrellaParameter> params;
};
