#include <vector> 
#include <string>
#include <mpi.h>
#include <lammps.h>
#include <input.h>
#include <atom.h>
#include <library.h>
#include "umbrella_step.h"

double *UmbrellaStep::positions_buffer = NULL;

UmbrellaStep::UmbrellaStep(LAMMPS_NS::LAMMPS *lmp, float probability, char* name) {
	this->lmp = lmp;
	this->probability = probability;
	strcpy (this->name, name);
	rand_min = 1;
	rand_max = 0;
	if (positions_buffer == NULL)
		positions_buffer = new double[3 * 1000];
};
UmbrellaStep::UmbrellaStep() {}

UmbrellaStep::~UmbrellaStep(){};

void UmbrellaStep::execute_init() {
	execute_block (lmp, step_init_block);
}
void UmbrellaStep::execute_step() {
	execute_block (lmp, take_step_block);
}

void UmbrellaStep::execute_accept() {
	execute_block (lmp, if_accept_block);
}

void UmbrellaStep::execute_reject() {
	execute_block (lmp, if_reject_block);
}

// static
void UmbrellaStep::execute_block (LAMMPS_NS::LAMMPS *lmp, std::vector<std::string> block) {
	for (int i = 0; i < block.size(); ++i) {
		if (strcmp (block[i].c_str(), "GET_ATOMS") == 0) {
			lmp->input->one ("# Copying positions to buffer...");
			lammps_gather_atoms(lmp, (char*)"x", 1, 3, positions_buffer);
		} else if (strcmp (block[i].c_str(), "PUT_ATOMS") == 0) {
			lmp->input->one ("# Scattering buffered positions...");
			lammps_scatter_atoms(lmp, (char*)"x", 1, 3, positions_buffer);
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
