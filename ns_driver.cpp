/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   www.cs.sandia.gov/~sjplimp/lammps.html
   Steve Plimpton, sjplimp@sandia.gov, Sandia National Laboratories

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under 
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

// c++_driver = simple example of how an umbrella program
//              can invoke LAMMPS as a library on some subset of procs
// Syntax: c++_driver P in.lammps
//         P = # of procs to run LAMMPS on
//             must be <= # of procs the driver code itself runs on
//         in.lammps = LAMMPS input script
// See README for compilation instructions

#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "mpi.h"
#include "lammps.h"         // these are LAMMPS include files
#include "input.h"
#include "atom.h"
#include "library.h"

#include <math.h>
#include <sys/time.h>
#include "parser.h"
#include "pmi.h"	//nid

using namespace LAMMPS_NS;

int64_t get_time();
int main(int narg, char **arg)
{
	// setup MPI and various communicators
	// driver runs on all procs in MPI_COMM_WORLD
	// comm_lammps only has 1st P procs (could be all or any subset)
	MPI_Init(&narg,&arg);
	parser *p;

	//TODO take care to minimize number of nodes working on a given window
	//MPI things ~~~~~~~~~~~~~~~~~~~~~~~~
	int me,nprocs;
	MPI_Comm_rank(MPI_COMM_WORLD,&me);
	MPI_Comm_size(MPI_COMM_WORLD,&nprocs);

	int nid;
	PMI_CNOS_Get_nid(me, &nid);

	p->verbose = me == 0;
	p = new parser (narg, arg);
	if (p->parse_error)
		return 1;

	const int window_index = me % p->num_windows;
	
	MPI_Comm local_comm;
	//TODO warn when nprocs % num_windows != 0, error when nprocs > num_windows
	MPI_Comm_split(MPI_COMM_WORLD, me % p->num_windows, me, &local_comm);
	int local_rank;
	MPI_Comm_rank (local_comm, &local_rank);

	//Setup logging for this window TODO window, not proc
	FILE *local_log;
	if (local_rank == 0) {
		char out_fname[100];
		sprintf (out_fname, "output/log_%d.txt", window_index);
		local_log = fopen (out_fname, "w");
		setbuf (local_log, NULL);//debug
	}

	int my_file_index;
	int file_indices[p->num_windows];//only the root process needs allocate this
	float targets[p->num_windows];
	if (me == 0)
	{
		FILE *in_file = fopen ("sn_window_init.txt", "r");
		for (int i = 0; i < p->num_windows; i++) {
			fscanf (in_file, "%d %f", &file_indices[i], &targets[i]);
		}
		fclose (in_file);
	}
	//TODO scatterv this, assumes num_windows procs, 1 proc/window
	MPI_Bcast (file_indices, p->num_windows, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast (targets, p->num_windows, MPI_FLOAT, 0, MPI_COMM_WORLD);


	if (local_rank == 0) {
		fprintf (stderr, "I am proc %d working on file index %d with target %f, running on node ID: %d.\n", me, file_indices[window_index], targets[window_index], nid);
		fprintf (local_log, "I am proc %d working on file index %d with target %f, running on node ID: %d.\n", me, file_indices[window_index], targets[window_index], nid);
	}

	
	int lammps = 1;//TODO clean this up
	LAMMPS *lmp;



	if (lammps == 1) {
		lmp = new LAMMPS(0,NULL,local_comm);
		int n;
		char line[1024];
		lmp->input->one("#LAMMPS umbrella sampler");
		lmp->input->one("#Anthony Frachioni 10-2012");
		lmp->input->one("#Hardcoded umbrella script for demo");
		lmp->input->one("units metal");
		//sprintf (line, "read_restart %s", p->start);
		sprintf (line, "read_data quench_init/input/input_sn.%d", file_indices[window_index]);
		lmp->input->one(line);
		lmp->input->one("reset_timestep 0");
		lmp->input->one("compute myBoop all boop 6 3.4");//TODO parse
		lmp->input->one("thermo 5");
		sprintf (line, "log output/log_%d.lammps", window_index);
		lmp->input->one(line);

		//Ar
		//lmp->input->one("pair_style lj/cut 5");
		//lmp->input->one("pair_coeff 1 1 0.0104 3.405");//TODO hardcoded

		//Sn
		lmp->input->one("pair_style meam");
		lmp->input->one("pair_coeff * * library.meam Sn NULL Sn");//TODO hardcoded

		lmp->input->one("run 0");
		double Q6 = *((double *) lammps_extract_compute(lmp, (char*)"myBoop", 0, 0));

		sprintf (line, "fix 1 all npt temp %f %f 0.1 iso 1.0 1.0 0.1", p->temperature, p->temperature);
		lmp->input->one(line);

		int natoms = lammps_get_natoms (lmp);//XXX uses 32 bit int
		if (local_rank == 0) {
			fprintf (local_log, "Initial Q6: %f\n", Q6);
			fprintf (local_log, "Number of atoms: %d\n", natoms);
			fprintf (local_log, "Spring/Boltzmann: %f\n", p->spring);
		}
		double *positions_buffer = new double[3 * natoms];
		double factor = -0.5 * p->spring / p->temperature;
		int umbrella_accept = 1;
		int seed;
		int accept_count = 0;
		int64_t last_step_end_time, step_time, lammps_start_time, lammps_split;
		last_step_end_time = get_time();
		double Q6_old = Q6;
		double d_Q, d_Q_old, bias_potential_new, bias_potential_old, log_boltz_factor;
		//Q6_old is last accepted Q6
		for (int i = 0; i < p->count; i++) {
			if (!umbrella_accept) {
				//Umbrella reject
				//TODO try copying lmp->atom->x to avoid global op.
				lammps_put_coords (lmp, positions_buffer);
				if (local_rank == 0) seed = rand() % 1000 + 1; //TODO why 1000?
				MPI_Bcast (&seed, 1, MPI_INT, 0, local_comm);
				sprintf(line, "velocity all create %f %d", p->temperature, seed);
				lmp->input->one (line);
			} else {
				//Umbrella accept
				accept_count++;
				lammps_get_coords (lmp, positions_buffer);
				Q6_old = Q6;
				d_Q_old = Q6_old - targets[window_index];
				bias_potential_old = d_Q_old * d_Q_old;
			}
			sprintf (line, "run %d", p->duration);
			lammps_start_time = get_time();
			lmp->input->one (line);
			lammps_split = get_time() - lammps_start_time;
			Q6 = *((double *) lammps_extract_compute(lmp,(char*)"myBoop", 0, 0));
			if (local_rank == 0) {
				d_Q = Q6 - targets[window_index];
				bias_potential_new = d_Q * d_Q;
				log_boltz_factor = factor*(bias_potential_new-bias_potential_old);
				umbrella_accept = log((double) rand() / RAND_MAX) < log_boltz_factor;
			}
			MPI_Bcast (&umbrella_accept, 1, MPI_INT, 0, local_comm);
			step_time = get_time() - last_step_end_time;
			if (local_rank == 0){
				fprintf (local_log, "[%d] Step index: %d\tQ6: %f\t (log)Boltzmann factor: %f\tAccept: %d\tAccept rate: %f\tTime: %5lld\tLAMMPS time: %5lld\tWall: %lld\n", window_index, i, Q6, log_boltz_factor, umbrella_accept, ((double) accept_count) / (i + 1), step_time, lammps_split, get_time());
				//fprintf (local_log, "\td_Q: %f\td_Q_old: %f\tbias: %f\tbias_old: %f\tfactor: %f\n", d_Q, d_Q_old, bias_potential_new, bias_potential_old, factor);
			}
			last_step_end_time = get_time();
		}
		delete [] positions_buffer;
		delete lmp;
	}
	fclose (local_log);
	MPI_Finalize();
}
int64_t get_time()
{
	struct timeval tv;
	gettimeofday(&tv, NULL);
	uint64_t ret = tv.tv_usec;
	ret /= 1000;
	ret += (tv.tv_sec * 1000);
	return ret;
}
