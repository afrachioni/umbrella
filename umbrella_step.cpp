#include <vector> 
#include <string>
#include <mpi.h>
#include <lammps.h>
#include <input.h>
#include <atom.h>
#include <library.h>
#include "umbrella_step.h"

double *UmbrellaStep::positions_buffer = NULL;
int UmbrellaStep::get_atoms_called = 0;
// The zeroeth step is always accepted
int UmbrellaStep::force_accept = 1;

UmbrellaStep::UmbrellaStep(LAMMPS_NS::LAMMPS *lmp, float probability, char* name, Global *global) {
	this->lmp = lmp;
	this->probability = probability;
	strcpy (this->name, name);
	rand_min = 1;
	rand_max = 0;
	this->global = global;
};
UmbrellaStep::UmbrellaStep() {}

UmbrellaStep::~UmbrellaStep(){};

void UmbrellaStep::execute_init() {
	execute_block (lmp, step_init_block, global);
}
void UmbrellaStep::execute_step() {
	execute_block (lmp, take_step_block, global);
}

void UmbrellaStep::execute_accept() {
	execute_block (lmp, if_accept_block, global);
}

void UmbrellaStep::execute_reject() {
	execute_block (lmp, if_reject_block, global);
}

// static
void UmbrellaStep::execute_block (LAMMPS_NS::LAMMPS *lmp, std::vector<std::string> block, Global *global) {
	for (int i = 0; i < block.size(); ++i) {
		if (strcmp (block[i].c_str(), "GET_ATOMS") == 0) {
			if (!get_atoms_called) {
				lmp->input->one ("run 0");
				positions_buffer = new double[3 * lmp->atom->natoms];
				get_atoms_called = 1;
			}
			lmp->input->one ("# Copying positions to buffer...");
			lammps_gather_atoms(lmp, (char*)"x", 1, 3, positions_buffer);
		} else if (strcmp (block[i].c_str(), "PUT_ATOMS") == 0) {
			if (!get_atoms_called)
				global->abort((char*) "put_positons was called "
						"before any call to get_positons.");
			lmp->input->one ("# Scattering buffered positions...");
			lammps_scatter_atoms(lmp, (char*)"x", 1, 3, positions_buffer);
		} else if (strcmp (block[i].c_str(), "FORCE_ACCEPT") == 0) {
			force_accept = 1;
		} else
			lmp->input->one (block[i].c_str());
	}
}

std::vector<std::string>* UmbrellaStep::get_take_step_block() {
	return &take_step_block;
}

std::vector<std::string>* UmbrellaStep::get_if_accept_block() {
	return &if_accept_block;
}

std::vector<std::string>* UmbrellaStep::get_if_reject_block() {
	return &if_reject_block;
}

std::vector<std::string>* UmbrellaStep::get_step_init_block() {
	return &step_init_block;
}
