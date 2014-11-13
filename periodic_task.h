#ifndef PERIODIC_TASK_H
#define PERIODIC_TASK_H

#include <mpi.h>
#include <lammps.h>

class PeriodicTask {
	public:
		PeriodicTask(LAMMPS_NS::LAMMPS *lmp, int period);

		std::vector<std::string>* get_task_block();

		int period;
		void execute_task(int step);

	private:
		LAMMPS_NS::LAMMPS *lmp;
		std::vector<std::string> task_block;
};
#endif
