#include<vector>
#include<string>
#include<mpi.h>
#include<lammps.h>

#include "periodic_task.h"
#include "umbrella_step.h"

PeriodicTask::PeriodicTask(LAMMPS_NS::LAMMPS *lmp, int period, Global *global){
		this->lmp = lmp;
		this->period = period;
		this->global = global;
}

void PeriodicTask::execute_task(int step) {
	if (step % period == 0)
		UmbrellaStep::execute_block (lmp, task_block, global);
}

std::vector<std::string>* PeriodicTask::get_task_block() {
	return &task_block;
}
