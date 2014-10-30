#include <vector> 
#include <string.h>
#include <mpi.h>
#include <lammps.h>
#include <input.h>
#include <atom.h>
#include <library.h>
#include "umbrella_step.h"

double *UmbrellaStep::positions_buffer = NULL;
int *UmbrellaStep::types_buffer = NULL;
double UmbrellaStep::xlo, UmbrellaStep::ylo, UmbrellaStep::zlo, \
UmbrellaStep::xhi, UmbrellaStep::yhi, UmbrellaStep::zhi;
int UmbrellaStep::get_atoms_called = 0;
int UmbrellaStep::get_types_called = 0;
char UmbrellaStep::line[100];
// The zeroeth step is always accepted
int UmbrellaStep::force_accept = 1;

UmbrellaStep::UmbrellaStep(LAMMPS_NS::LAMMPS *lmp, float probability, char* name, Global *global) {
	is_barostat = 0;
	this->lmp = lmp;
	this->probability = probability;
	strcpy (this->name, name);
	rand_min = 1;
	rand_max = 0;
	this->global = global;

	this->logger = NULL;
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
	for (unsigned i = 0; i < block.size(); ++i) {
		if (strcmp (block[i].c_str(), "GET_ATOMS") == 0) {
			if (!get_atoms_called) {
				lmp->input->one ("run 0");
				positions_buffer = new double[3 * lmp->atom->natoms];
				get_atoms_called = 1;
			}
			lmp->input->one ("# Copying positions to buffer...");
			lammps_gather_atoms(lmp, (char*)"x", 1, 3, positions_buffer);
			// Perhaps it's safe to copy, modify internal box bounds array
			// Cursory check: it probably is.
			xlo = *((double *) lammps_extract_global(lmp, (char*)"boxxlo"));
			ylo = *((double *) lammps_extract_global(lmp, (char*)"boxylo"));
			zlo = *((double *) lammps_extract_global(lmp, (char*)"boxzlo"));
			xhi = *((double *) lammps_extract_global(lmp, (char*)"boxxhi"));
			yhi = *((double *) lammps_extract_global(lmp, (char*)"boxyhi"));
			zhi = *((double *) lammps_extract_global(lmp, (char*)"boxzhi"));

			sprintf (line, "# Positions retrieved. 1: %f\t%f\t%f",
					positions_buffer[0], positions_buffer[1],
					positions_buffer[2]);
			lmp->input->one (line);

			sprintf (line, "# Positions retrieved. 2: %f\t%f\t%f",
					positions_buffer[3], positions_buffer[4],
					positions_buffer[5]);
			lmp->input->one (line);

		} else if (strcmp (block[i].c_str(), "PUT_ATOMS") == 0) {
			if (!get_atoms_called)
				global->abort((char*) "put_positons was called "
						"before any call to get_positons.");
			lmp->input->one ("# Scattering buffered positions...");

			sprintf (line, "change_box all x final %f %f units box", xlo, xhi);
			lmp->input->one (line);
			sprintf (line, "change_box all y final %f %f units box", ylo, yhi);
			lmp->input->one (line);
			sprintf (line, "change_box all z final %f %f units box", zlo, zhi);
			lmp->input->one (line);

			int64_t time = Logger::get_time();
			lammps_scatter_atoms(lmp, (char*)"x", 1, 3, positions_buffer);
			if (Global::get_instance()->get_global_rank() == 0) {
				fprintf (stdout, "time: %" PRId64 "\n", Logger::get_time());
				fprintf (stdout, "scatter time: %" PRId64 "\n", Logger::get_time() - time);
			}

			sprintf (line, "# Positions scattered. 1: %f\t%f\t%f",
					positions_buffer[0], positions_buffer[1],
					positions_buffer[2]);
			lmp->input->one (line);

			sprintf (line, "# Positions scattered. 2: %f\t%f\t%f",
					positions_buffer[3], positions_buffer[4],
					positions_buffer[5]);
			lmp->input->one (line);


		} else if (strcmp (block[i].c_str(), "GET_TYPES") == 0) {
			if (!get_types_called) {
				lmp->input->one ("run 0"); //XXX why?
				types_buffer = new int[lmp->atom->natoms];
				get_types_called = 1;
			}
			lmp->input->one ("# Copying types to buffer...");
			lammps_gather_atoms(lmp, (char*)"type", 0, 1, types_buffer);
		} else if (strcmp (block[i].c_str(), "PUT_TYPES") == 0) {
			if (!get_types_called)
				global->abort((char*) "put_types was called "
						"before any call to get_positions.");
			lmp->input->one("# Scatternig buffered types...");;
			lammps_scatter_atoms(lmp, (char*)"type", 0, 1, types_buffer);
		} else if (strcmp (block[i].c_str(), "FORCE_ACCEPT") == 0) {
			force_accept = 1;
		//} else if (strcmp (block[i].c_str(), "DO_STEP") == 0) {
			// execute step
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

void UmbrellaStep::set_logger_debug (Logger* logger) {
	this->logger = logger;
}
