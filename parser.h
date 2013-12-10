#include <map>
#include <mpi.h>
#include <lammps.h>

class Parser {
	public:
		Parser(const char *fname, LAMMPS_NS::LAMMPS *lmp);
		~Parser();
		void parse();
		void print();
	private:
		char *fname;
		LAMMPS_NS::LAMMPS *lmp;
		std::vector<std::string> init_block;
		std::map<std::string, UmbrellaStep> steps_map;
		std::vector<UmbrellaParameter> params;
};
