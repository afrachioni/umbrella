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

#define VERSION "13.02.02.1"

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
	
	if (me == 0) {
		fprintf (stderr, "                       -------------------------------------------------\n");
		fprintf (stderr, "                       |  LAMMPS umbrella sampler, version %s  |\n", VERSION);
		fprintf (stderr, "                       |           Anthony Frachioni,  1-2013          |\n");
		fprintf (stderr, "                       |            afrachi1@binghamton.edu            |\n");
		fprintf (stderr, "                       -------------------------------------------------\n");
	}

	int nid;
	PMI_CNOS_Get_nid(me, &nid);

	p->verbose = me == 0;
	p = new parser (narg, arg);
	if (p->parse_error)
		return 1;

	//read init file once to determine number of windows
	int num_windows = 0;
	int num_lines = 0;
	char c = 'p';
	char d;
	c = fgetc (p->init);//these three lines only care if first line is commented
	ungetc (c, p->init);
	if (c == '#') --num_windows;
	while (c != EOF) {
		c = fgetc (p->init);
		if (c == '\n') {
			++num_lines;
			d = fgetc (p->init);
			ungetc (d, p->init);
			if (d != '#')
				++num_windows;
			d = 'p';
		}
	}
	rewind (p->init);
	//fprintf (stderr, "It looks like there are %d lines and %d windows in the init file.\n", num_lines, num_windows);

	if (me == 0) {
		if (nprocs < num_windows)
			fprintf (stderr, "WARNING: More windows have been defined (%d) than exist available processors (%d).  The first %d windows will execute.\n", num_windows, nprocs, nprocs);
		else if (nprocs % num_windows != 0)
			fprintf (stderr, "WARNING: The number of available processors (%d) is not evenly divisible by the number of defined windows (%d).\n"
					"         Windows 0 to %d will use %d processors, windows %d to %d will use %d processors.\n", nprocs, num_windows, nprocs % num_windows - 1, nprocs / num_windows + 1, nprocs % num_windows, num_windows - 1, nprocs / num_windows);
		else if (p->verbose)
			fprintf (stderr, "There are %d available processors and %d defined windows; each window will use %d processors.\n", nprocs, num_windows, nprocs / num_windows);
	}

	const int window_index = me % num_windows;
	
	//Split COMM_WORLD into communicators for each window
	MPI_Comm local_comm;
	//TODO warn when nprocs % num_windows != 0, error when nprocs > num_windows
	MPI_Comm_split(MPI_COMM_WORLD, me % num_windows, me, &local_comm);
	int local_rank;
	MPI_Comm_rank (local_comm, &local_rank);


	//Setup logging for this window
	FILE *local_log;
	if (local_rank == 0) {
		char out_fname[100];
		sprintf (out_fname, "output/log_%d.txt", window_index);
		local_log = fopen (out_fname, "w");
		setbuf (local_log, NULL);//debug
	}

	//Read init file again, this time store data
	char file_paths[num_windows][200];//only root process
	float targets[num_windows];
	if (me == 0)
	{
		int j = 0;
		for (int i = 0; i < num_lines; i++) {
			c = fgetc (p->init);
			if (c == '#') {
				fscanf (p->init, "%*99[^\n]\n");
				continue;
			}
			ungetc (c, p->init);
			fscanf (p->init, "%s %f\n", &file_paths[j], &targets[j]);
			//fprintf (stderr, "Root has read j = %2d, c = %c, file_paths[j] = %s, targets[j] = %f\n", j, c, file_paths[j], targets[j]);
			++j;
		}
		fclose (p->init);
		if (j != num_windows) fprintf (stderr, "Parse error! j: %d, num_windows: %d\n", j, num_windows);


		//write meteadada file for WHAM
		FILE *meta_file = fopen ("wham_metadata.txt", "w");
		fprintf (meta_file, "# Metadata file for Alan Grossfield's WHAM\n");
		fprintf (meta_file, "# Generated by Anthony Frachioni's LAMMPS umbrella sampler, version %s\n", VERSION);
		fprintf (meta_file, "#\n");
		for (int i = 0; i < num_windows; ++i)
			fprintf (meta_file, "output/log_%d.txt\t%f\t%f\n", i, targets[i], p->spring * WHAM_BOLTZMANN);
		fclose (meta_file);
	}
	//TODO scatter to local roots, than bcast internally
	MPI_Bcast (file_paths, num_windows * 200, MPI_CHAR, 0, MPI_COMM_WORLD);
	MPI_Bcast (targets, num_windows, MPI_FLOAT, 0, MPI_COMM_WORLD);


	//Setup LAMMPS instance with initial configuration, proper units, Q compute, etc.
	LAMMPS *lmp = new LAMMPS(0,NULL,local_comm);
	int n;
	char line[100];
	lmp->input->one("#LAMMPS umbrella sampler");
	lmp->input->one("#Anthony Frachioni 10-2012");
	lmp->input->one("#Hardcoded umbrella script for demo");
	lmp->input->one("units metal");
	sprintf (line, "read_data %s", file_paths[window_index]);
	lmp->input->one(line);
	lmp->input->one("reset_timestep 0");
	sprintf (line, "compute myBoop all boop %d %f", p->l, p->cutoff);
	lmp->input->one(line);
	lmp->input->one("thermo 5");
	sprintf (line, "log output/log_%d.lammps", window_index);
	lmp->input->one(line);
	lmp->input->one("pair_style meam");
	lmp->input->one("pair_coeff * * library.meam Sn NULL Sn");//TODO hardcoded potential
	lmp->input->one("run 0");
	sprintf (line, "fix 1 all npt temp %f %f 0.1 iso 1.0 1.0 0.1", p->temperature, p->temperature);
	lmp->input->one(line);
	double Q6 = *((double *) lammps_extract_compute(lmp, (char*)"myBoop", 0, 0));
	int natoms = lammps_get_natoms (lmp);//XXX library uses 32 bit int for some reason

	//Write header information to local logfile
	if (local_rank == 0) {
		fprintf (local_log, "# Logfile from Anthony Frachioni's LAMMPS umbella sampler\n");
		fprintf (local_log, "# Should be directly readable by Grossfield's WHAM code\n");
		fprintf (local_log, "# Generated using version: %s\n", VERSION);
		fprintf (local_log, "# I am proc %d working on window index %d, running on node ID: %d.\n", me, window_index, nid);
		fprintf (local_log, "# Using initial configuration: %s\n", file_paths[window_index]);
		fprintf (local_log, "# Target Q6: %f\n", targets[window_index]);//TODO change "Q6" to something, merge these prints with those below
		fprintf (local_log, "# Initial Q6: %f\n", Q6);
		fprintf (local_log, "# Number of atoms: %d\n", natoms);
		fprintf (local_log, "# Spring/Boltzmann: %f\n", p->spring);
		fprintf (local_log, "#\n");
		//fprintf (local_log, "#%-14s\t%-15s\t%-15s\t%-8s\t%-15s\t%-15s\n", "Step", "Last", "Current", "Accept", "LAMMPS_time", "Wall");
		fprintf (local_log, "#%-8s%-15s%-15s%-8s%-15s%-15s\n", "Step", "Last", "Current", "Accept", "LAMMPS_time", "Wall");
		//formerly 14 here   ^^
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
			//TODO try memcopying lmp->atom->x to avoid global op.
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
			//fprintf (local_log, "%-15d\t%-15f\t%-15f\t%-8d\t%-15lld\t%-15lld\n", i, Q6_old, Q6, umbrella_accept, lammps_split, get_time());
			fprintf (local_log, "%-9d%-15f%-15f%-8d%-15lld%-15lld\n", i, Q6_old, Q6, umbrella_accept, lammps_split, get_time());
			//formerly 15 here   ^^
			//fprintf (local_log, "%d\t%f\n", i, Q6_old);
		}
		last_step_end_time = get_time();
	}
	delete [] positions_buffer;
	delete lmp;
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
