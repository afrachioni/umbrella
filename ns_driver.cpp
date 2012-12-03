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

using namespace LAMMPS_NS;

int64_t get_time();
int main(int narg, char **arg)
{

	// setup MPI and various communicators
	// driver runs on all procs in MPI_COMM_WORLD
	// comm_lammps only has 1st P procs (could be all or any subset)
	MPI_Init(&narg,&arg);
	parser *p;

	//MPI things ~~~~~~~~~~~~~~~~~~~~~~~~
	int me,nprocs;
	MPI_Comm_rank(MPI_COMM_WORLD,&me);
	MPI_Comm_size(MPI_COMM_WORLD,&nprocs);

	int nprocs_lammps = nprocs;//TODO parse or ensemble or something
	if (nprocs_lammps > nprocs) {
		if (me == 0)
			printf("ERROR: LAMMPS cannot use more procs than available\n");
		MPI_Abort(MPI_COMM_WORLD,1);
	}
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	p->verbose = me == 0;
	p = new parser (narg, arg);
	if (p->parse_error)
		return 1;
	
	int lammps;
	if (me < nprocs_lammps) lammps = 1;
	else lammps = MPI_UNDEFINED;
	MPI_Comm comm_lammps;
	MPI_Comm_split(MPI_COMM_WORLD,lammps,0,&comm_lammps);
	LAMMPS *lmp;



	if (lammps == 1) {
		lmp = new LAMMPS(0,NULL,comm_lammps);
		int n;
		char line[1024];
		lmp->input->one("#LAMMPS umbrella sampler");
		lmp->input->one("#Anthony Frachioni 10-2012");
		lmp->input->one("#Hardcoded umbrella script for demo");
		lmp->input->one("units metal");
		//sprintf (line, "read_restart %s", p->start);
		sprintf (line, "read_data %s", p->start);
		lmp->input->one(line);
		lmp->input->one("reset_timestep 0");
		lmp->input->one("compute myBoop all boop 6 3.4");//TODO parse
		lmp->input->one("thermo 5");

		//Ar
		//lmp->input->one("pair_style lj/cut 5");
		//lmp->input->one("pair_coeff 1 1 0.0104 3.405");//TODO hardcoded

		//Sn
		lmp->input->one("pair_style meam");
		//XXX BSn?
		lmp->input->one("pair_coeff * * library.meam Sn NULL Sn");//TODO hardcoded

		lmp->input->one("run 0");
		double Q6 = *((double *) lammps_extract_compute(lmp, (char*)"myBoop", 0, 0));

		sprintf (line, "fix 1 all npt temp %f %f 0.1 iso 1.0 1.0 0.1", p->temperature, p->temperature);
		lmp->input->one(line);

		int natoms = lammps_get_natoms (lmp);//XXX uses 32 bit int
		if (me == 0) {
			fprintf (stderr, "Initial Q6: %f\n", Q6);
			fprintf (stderr, "Number of atoms: %d\n", natoms);
		}
		double *positions_buffer = new double[3 * natoms];
		double factor = -0.5 * p->spring / p->temperature;
		int umbrella_accept = 1;
		int seed;
		int accept_count = 0;
		int64_t last_step_end_time, step_time;
		int64_t lammps_start_time, lammps_split;
		last_step_end_time = get_time();
		double Q6_old = Q6;
		double d_Q, d_Q_old, bias_potential_new, bias_potential_old;
		double boltz_factor;
		for (int i = 0; i < p->count; i++) {
			if (!umbrella_accept) {
				//TODO try copying lmp->atom->x to avoid global op.
				lammps_get_coords (lmp, positions_buffer);
				//TODO why 1000?
				if (me == 0) seed = rand() % 1000;
				MPI_Bcast (&seed, 1, MPI_INT, 0, comm_lammps);
				sprintf(line, "velocity all create %f %d", p->temperature, seed);
				lmp->input->one (line);
			} else {
				accept_count++;
				lammps_put_coords (lmp, positions_buffer);
			}
			sprintf (line, "run %d", p->duration);
			lammps_start_time = get_time();
			lmp->input->one (line);
			lammps_split = get_time() - lammps_start_time;
			Q6 = *((double *) lammps_extract_compute(lmp,(char*)"myBoop", 0, 0));
			if (me == 0) {
				d_Q_old = Q6_old - p->target_min;
				d_Q = Q6 - p->target_min;
				bias_potential_old = d_Q_old * d_Q_old;
				bias_potential_new = d_Q * d_Q;
				boltz_factor = exp(factor*(bias_potential_old-bias_potential_new));
				umbrella_accept = (((double) rand()) / RAND_MAX) < boltz_factor;
			}
			MPI_Bcast (&umbrella_accept, 1, MPI_INT, 0, comm_lammps);
			Q6_old = Q6;
			step_time = get_time() - last_step_end_time;
			if (me == 0)
				fprintf (stderr, "Step index: %d\tQ6: %f\tBoltzmann factor: %f\tAccept: %d\tAccept rate: %f\tTime: %5lld\tLAMMPS time: %5lld\tOverhead: %f\n", i, Q6, boltz_factor, umbrella_accept, ((double) accept_count) / (i + 1), step_time, lammps_split, ((double ) (step_time - lammps_split)) / step_time);
			last_step_end_time = get_time();
		}
		delete [] positions_buffer;
		delete lmp;
	}
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
