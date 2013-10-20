/* ----------------------------------------------------------------------
                     LAMMPS umbrella sampling driver
   Links against a modified LAMMPS to execute umbrella sampling
   Uses custom "compute boop" to calculate order parameter in situ
   Built to examine nucleation in tin

   Copyright 2013 Anthony Frachioni - All Rights Reserved
                           afrachi1@binghamton.edu
-------------------------------------------------------------------------
   This program is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program.  If not, see <http://www.gnu.org/licenses/>.
------------------------------------------------------------------------- */

#define VERSION "13.10.7.0"
#define DEBUG 1
#define MAX_FNAME_LENGTH 500
#define DUMP_EVERY_STEPS 100

//Command line options parser
#include "parser.h"


//LAMMPS include files
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <mpi.h>
#include <lammps.h>
#include <input.h>
#include <atom.h>
#include <library.h>
#include "management.cpp"

//system include files
#include <math.h>
#include <pthread.h>
#include <sys/time.h>
#include <sys/stat.h> //mkdir needs this
#include <sys/types.h> //and this too

//print node ID on Kraken
#ifdef KRAKEN
#include <pmi.h>
#endif

#define printmsg(...) if (me == 0) fprintf(stdout, __VA_ARGS__);
#define debugmsg(...) if (DEBUG && me == 0) fprintf(stdout, __VA_ARGS__);

using namespace LAMMPS_NS;

int64_t get_time();
int main(int narg, char **arg)
{
		MPI_Init(&narg,&arg);
		parser *p;
		int me,nprocs;
		MPI_Comm_rank(MPI_COMM_WORLD,&me);
		MPI_Comm_size(MPI_COMM_WORLD,&nprocs);
		setbuf (stdout, NULL);
		setbuf (stderr, NULL);

		printmsg ("                       "
				"-------------------------------------------------\n");
		printmsg ("                       "
				"|  LAMMPS umbrella sampler, version %s  |\n", VERSION);
		printmsg ("                       "
				"|           Anthony Frachioni,  1-2013          |\n");
		printmsg ("                       "
				"|            afrachi1@binghamton.edu            |\n");
		printmsg ("                       "
				"-------------------------------------------------\n");

		//get node ID if on Kraken
		int nid  = -1;
#ifdef KRAKEN
		PMI_CNOS_Get_nid(me, &nid);
#endif

		//parse command line arguments
		p->verbose = me == 0;//turn off for DEBUG?
		p = new parser (narg, arg);
		if (p->parse_error) {
			printmsg ("Parse errors present, exiting...\n");
			MPI_Finalize ();
		}

		//read init file once to determine number of windows
		//FIXME only root process should do this
		int num_windows = 0;
		int num_lines = 0;
		char c = 'p';
		char d;
		//if global rank = 0
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
		//bcast num_windows to everybody.
		//now why was that so hard?
		printmsg ("\n");

		if (nprocs < num_windows) {
			printmsg ("WARNING: More windows have been defined (%d) than "
					"exist available processors (%d).  The first %d windows "
					"will execute.\n", num_windows, nprocs, nprocs);
		} else if (nprocs % num_windows != 0) {
			printmsg ("WARNING: The number of available processors (%d) is "
					"not evenly divisible by the number of defined windows "
					"(%d).\n         Windows 0 to %d will use %d "
					"processors, windows %d to %d will use %d "
					"processors.\n", nprocs, num_windows, \
					nprocs % num_windows - 1, nprocs / num_windows + 1, \
					nprocs % num_windows, num_windows - 1, \
					nprocs / num_windows);
		} else if (p->verbose) {
			printmsg ("There are %d available processors and %d defined "
					"windows; each window will use %d processors.\n", \
					nprocs, num_windows, nprocs / num_windows);
		}

		int num_active_windows;
		if (nprocs < num_windows)
			num_active_windows = nprocs;
		else
			num_active_windows = num_windows;

		//Split COMM_WORLD into communicators for each window
		MPI_Comm local_comm;
		MPI_Comm_split(MPI_COMM_WORLD, me % num_active_windows, me, &local_comm);
		int local_rank;
		MPI_Comm_rank (local_comm, &local_rank);
		const int window_index = me % num_active_windows;
		int local_roots_send[num_active_windows];
		int local_roots[num_active_windows];
		for (int i = 0; i < num_active_windows; ++i)
			local_roots_send[i] = 0;
		if (local_rank == 0)
			local_roots_send[window_index] = me;
			MPI_Allreduce (local_roots_send, local_roots, num_active_windows, \
					MPI_INT, MPI_SUM, MPI_COMM_WORLD);
		//if (me == 0)
			//for (int i = 0; i < num_active_windows; ++i)
				//fprintf(stdout, "Root of window %2d is proc %2d.\n", i, local_roots[i]);

		MPI_Group world_group, roots_group;
		MPI_Comm_group (MPI_COMM_WORLD, &world_group);
		//fprintf (stdout, "rank: %2d, num_active_windows: %d\n", me, num_active_windows);
		//fprintf (stdout, "rank: %2d, local_roots[0]: %d\n", me, local_roots[0]);
		//fprintf (stdout, "rank: %2d, local_roots[1]: %d\n", me, local_roots[1]);
		//fprintf (stdout, "rank: %2d, local_roots[2]: %d\n", me, local_roots[2]);
		//fprintf (stdout, "rank: %2d, local_roots[3]: %d\n", me, local_roots[3]);
		MPI_Group_incl (world_group, num_active_windows, local_roots, \
				&roots_group);
		MPI_Comm roots_comm;
		MPI_Comm_create (MPI_COMM_WORLD, roots_group, &roots_comm);

		debugmsg ("Opening window logfiles...\n");
		//Setup logging for this window
		FILE *local_log;
		if (local_rank == 0) {
			char out_fname[100];
			mkdir ("logs", 0755);
			sprintf (out_fname, "logs/log_%d.txt", window_index);
			local_log = fopen (out_fname, "w");
			setbuf (local_log, NULL);//debug
		}

		debugmsg ("Reading window initialization data...\n");
		//Read init file again, this time store data
		//TODO move this block just after rewind (p->init); maybe
		char file_paths[num_active_windows][MAX_FNAME_LENGTH];
		float targets[num_active_windows];
		float springs[num_active_windows];
		if (me == 0) {
			int j = 0;
			for (int i = 0; i < num_lines && j < num_active_windows; i++) {
				c = fgetc (p->init);
				if (c == '#') {
					fscanf (p->init, "%*99[^\n]\n");
					continue;
				}
				ungetc (c, p->init);
				fscanf (p->init, "%s %f %f\n", \
						&file_paths[j], &targets[j], &springs[j]);
				++j;
			}
			fclose (p->init);
			//if (j != num_windows) 
				//printmsg ("Parse error! j: %d, num_windows: %d\n", \
						//j, num_windows);

			debugmsg ("Writing WHAM metadata file...\n");
			//write meteadada file for WHAM
			FILE *meta_file = fopen ("wham_metadata.txt", "w");
			fprintf(meta_file, "# Metadata file for Alan Grossfield's WHAM\n");
			fprintf (meta_file, "# Generated by Anthony Frachioni's LAMMPS "
					"umbrella sampler, version %s\n", VERSION);
			fprintf (meta_file, "#\n");
			for (int i = 0; i < num_active_windows; ++i)
				fprintf (meta_file, "logs/log_%d.txt\t%f\t%f\n", \
						i, targets[i], springs[window_index] * WHAM_BOLTZMANN);
			fclose (meta_file);
		}

		debugmsg ("Broadcasting initialization data...\n");
		//TODO scatter to local roots, than bcast internally
		MPI_Bcast (file_paths, num_active_windows * MAX_FNAME_LENGTH, MPI_CHAR, 0, \
				MPI_COMM_WORLD);
		MPI_Bcast (targets, num_active_windows, MPI_FLOAT, 0, MPI_COMM_WORLD);
		MPI_Bcast (springs, num_active_windows, MPI_FLOAT, 0, MPI_COMM_WORLD);
		
		//MPI_Scatter (springs, 1, MPI_FLOAT,
				//&local_spring, 1, MPI_FLOAT, 0, MPI_COMM_WORLD);
		//float *local_spring = new float (springs[window_index]);
		float spring_init = springs[window_index];

		//if (local_rank == 0)
			//fprintf (stderr, "window_index: %d\tspring_init: %f\n", window_index, spring_init);
		///////////////////////////////////

		char line[100];
		const int dump_freq = p->duration * DUMP_EVERY_STEPS;//TODO fix this in context of variable duration


		debugmsg ("Creating LAMMPSes...\n");
		sprintf (line, "logs/log_%d.screen", window_index);
		//Setup LAMMPS instance with initial conditions and settings
		//char *args[] = {(char*)"foo", (char*)"-screen", \
		line, (char*)"-log", (char*)"none"};
		char *args[] = {(char*)"foo", (char*)"-screen", (char*)"none", \
			(char*)"-log", (char*)"none"};
		LAMMPS *lmp = new LAMMPS(5,args,local_comm);
		//LAMMPS *lmp = new LAMMPS(0,NULL,local_comm);//for all stdout
#if DEBUG	
		sprintf (line, "log logs/log_%d.lammps", window_index);
		lmp->input->one(line);
#endif
		//It is necessary to create an atom map for usage of lammps_scatter_atoms.
		//Use hash map, assuming large number of atoms per processor
		lmp->input->one("atom_modify map hash");
		lmp->input->one("units metal");
		MPI_Barrier (MPI_COMM_WORLD);
		debugmsg ("About to read_data...\n");
		sprintf (line, "read_data %s", file_paths[window_index]);
		lmp->input->one(line);
		MPI_Barrier (MPI_COMM_WORLD);
		debugmsg ("Finished with read_data...\n");
		lmp->input->one("reset_timestep 0");
		sprintf (line, "compute Q6 all boop %d %f", p->l, p->cutoff);
		debugmsg ("About to call: %s\n", line);
		lmp->input->one(line);
		sprintf (line, "compute q6 all boop/atom %d %f", p->l, p->cutoff);
		debugmsg ("About to call: %s\n", line);
		lmp->input->one(line);
		debugmsg ("Finished with boop\n");
		lmp->input->one("pair_style meam");
		debugmsg ("About to call pair_coeff\n");
		lmp->input->one("pair_coeff * * library.meam Sn NULL Sn");//TODO hardcoded potential
		//lmp->input->one("pair_coeff * * library.meam Sn Ag NULL Sn Ag");
		lmp->input->one("thermo 5");
		debugmsg ("Defining thermo_style...\n");
		lmp->input->one("thermo_style custom step temp vol pe press c_Q6");
		debugmsg ("About to run 0...\n");
		lmp->input->one("run 0");
		debugmsg ("About to define fix...\n");
		//sprintf (line, "fix 1 all npt temp %f %f 0.1 iso 1.0 1.0 0.1", p->temperature, p->temperature);
		sprintf (line, "fix 1 all nve");
		lmp->input->one(line);
		debugmsg ("About to create dump and logs directories...\n");
		mkdir ("dump", 0755); mkdir ("logs", 0755);
		debugmsg ("About to define dump...\n");
		sprintf (line, "dump 1 all custom %d dump/dump_%d_*.txt id x y z c_q6", \
				dump_freq, window_index);
		//sprintf (line, "dump 1 all atom %d dump/dump_%d_*.txt", \
				//dump_freq, window_index);
		lmp->input->one(line);

		debugmsg ("Defining LAMMPS volume variable...\n");
		lmp->input->one("variable V equal vol"); //TODO determine why these are necessary
		lmp->input->one("variable U equal pe");
		debugmsg ("About to extract Q6...\n");
		double Q6 = *((double *) lammps_extract_compute(lmp, (char*)"Q6", 0, 0));
		double V = *((double *) lammps_extract_variable(lmp, (char*)"V", (char*)"all"));
		double U = *((double *) lammps_extract_variable(lmp, (char*)"U", (char*)"all"));
		debugmsg ("Initial volume: %f\n", V);
		debugmsg ("Initial potential energy: %f\n", U);
		int natoms = static_cast<int> (lmp->atom->natoms);

		debugmsg ("Writing header to window logs...\n");
		//Write header information to local logfile
		if (local_rank == 0) {
			fprintf (local_log, "# Logfile from Anthony Frachioni's LAMMPS umbella sampler\n");
			fprintf (local_log, "# Should be directly readable by Grossfield's WHAM code\n");
			fprintf (local_log, "# Generated using version: %s\n", VERSION);
			fprintf (local_log, "# I am proc %d working on window index %d, running on node ID: %d.\n", me, window_index, nid);
			fprintf (local_log, "# Using initial configuration: %s\n", file_paths[window_index]);
			fprintf (local_log, "# Target Q6: %f\n", targets[window_index]);//TODO change "Q6" to something, merge these prints with those below
			fprintf (local_log, "# Initial Q6: %f\n", Q6);
			fprintf (local_log, "# Initial potential energy / eV: %f\n", U);
			fprintf (local_log, "# Initial volume / A^3: %f\n", V);
			fprintf (local_log, "# Number of atoms: %d\n", natoms);
			fprintf (local_log, "# Spring/Boltzmann: %f\n", spring_init);//springs[window_index]);
			fprintf (local_log, "#\n");
			fprintf (local_log, "#%-11s%-15s%-15s%-8s%-15s%-15s%-8s%-15s\n", \
					"Step", "Last", "Current", "Accept", "LAMMPS_time", \
					"Wall", "Duration", "Spring");
		}

		debugmsg ("Writing WHAM command...\n");
		double target_max = 0;
		double target_min = INFINITY;
		for (int i = 0; i < num_active_windows; ++i) {
			if (target_max < targets[i])
				target_max = targets[i];
			if (target_min > targets[i])
				target_min = targets[i];
		}
		double target_lbound = target_min - 0.5 * (target_max - target_min) / num_active_windows;
		double target_rbound = target_max + 0.5 * (target_max - target_min) / num_active_windows;
		printmsg ("\nIt might be useful to analyze these data by running WHAM "
				"something like this:\n");
		sprintf (line, "wham %f %f %d 0.1 %f 0 wham_metadata.txt wham_free.txt", \
				target_lbound, target_rbound, 10 * num_active_windows, p->temperature);
		printmsg ("%s\n", line);
		printmsg ("which can also be accomplished by simply "
				"calling \'. wham_cmd\'.\n\n");
		FILE *wham_cmd_file = fopen ("wham_cmd.txt", "w");
		fprintf (wham_cmd_file, line);
		fclose (wham_cmd_file);


		//Umbrella definitions
		double *positions_buffer = new double[3 * natoms];
		int accept = 1;
		int md; // 1: MD move, 0: VMC move
		int seed;
		int accept_count = 0;
		int vmc_accept_count = 0;
		int64_t last_step_end_time, step_time, lammps_start_time, lammps_split;
		last_step_end_time = get_time();
		double Q6_old = Q6;
		double d_Q, d_Q_old, bias_potential_new, bias_potential_old;

		double V_old, U_old, d_V, d_U, dx, rec;
		const double P = 1.01325; //bar XXX this controls target pressure, should get passed in
		//const double pv_factor = 1.458397387e-5;// kcal / (atm A^3 mol) yes its wierd
		const double pv_factor = 6.0221413e-5;//kJ / (bar A^3 mol)
		const double du_factor = 9.648533632e1;//kJ / (eV mol)
		double log_boltz_factor;
		const int64_t loop_start_time = get_time();

		int duration_init = p->duration;
		int current_duration = duration_init;

		lammps_gather_atoms(lmp,(char*)"x",1,3,positions_buffer);
		Q6_old = Q6;
		d_Q_old = Q6_old - targets[window_index];
		bias_potential_old = d_Q_old * d_Q_old;
		U_old = U;
		V_old = V;

		pthread_mutex_t mpi_mutex = PTHREAD_MUTEX_INITIALIZER;

		const int count = 4;
		int lengths[count] = {1, 1, 1, 1};
		MPI_Aint offsets[count] = {0, sizeof(int), 2*sizeof(int), 3*sizeof(int)};
		MPI_Datatype types[count] = {MPI_INT, MPI_INT, MPI_INT, MPI_FLOAT};
		MPI_Datatype message_type;
		MPI_Type_struct (count, lengths, offsets, types, &message_type);
		MPI_Type_commit (&message_type);

		if (me == 0) {
			struct parameter_pointers management_data;
			management_data.window_index = window_index;
			management_data.num_active_windows = num_active_windows; 
			management_data.roots_comm = roots_comm;
			management_data.mpi_mutex_ptr = &mpi_mutex;
			management_data.message_type = message_type;

			pthread_t manager;
			if (pthread_create( &manager, NULL, &manage, &management_data))
				printmsg ("Could not create management thread on root!");
		}

		/*  XXX What's this?
			const int update_count = 4;
			int update_lengths[update_count] = {1, 1, 1, 1};
			MPI_Aint update_offsets[update_count] = {0, sizeof(int), 2*sizeof(int), 2*sizeof(int) + sizeof(double)};
			MPI_Datatype update_types[update_count] = {MPI_INT, MPI_INT, MPI_DOUBLE, MPI_DOUBLE};
			MPI_Datatype update_message_type;
			MPI_Type_struct (update_count, update_lengths, update_offsets, update_types, &update_message_type);
			MPI_Type_commit (&update_message_type);
		 */

		//mkfifo ("logs/feed_pipe", S_IWUSR | S_IRUSR);
		//int fd = open ("logs/feed_pipe", O_WRONLY | O_NONBLOCK);



		float current_spring = spring_init;
		MPI_Request req;
		int test_msg;
		struct management_message msg;
		int recv_complete = 0;
		if (local_rank == 0) {
			pthread_mutex_lock (&mpi_mutex);
			MPI_Irecv (&msg, 1, message_type, 0, TAG, roots_comm, &req);
			pthread_mutex_unlock (&mpi_mutex);
		}
		MPI_Barrier (MPI_COMM_WORLD);
		printmsg ("Samples away!\n\n");
		//Q6_old is most recently accepted Q6
		for (int i = 0; i < p->count; i++) {
			pthread_mutex_lock (&mpi_mutex);
			if (local_rank == 0) {
				MPI_Test (&req, &recv_complete, MPI_STATUS_IGNORE);
				if (recv_complete) {
					fprintf (stderr, "Window %d has recieved a message!\n", window_index);
					fprintf (stderr, "\tspring_msg: %d\n", msg.spring_message);
					fprintf (stderr, "\tduration_msg: %d\n", msg.duration_message);
					fprintf (stderr, "\tnew_duration: %d\n", msg.new_duration);
					fprintf (stderr, "\tnew_spring: %f\n", msg.new_spring);
					if (msg.spring_message)
						current_spring = msg.new_spring;
					if (msg.duration_message)
						current_duration = msg.new_duration;
					recv_complete = 0;
					MPI_Irecv (&msg, 1, message_type, 0, TAG, roots_comm, &req);
				}
			}
			//bcast to window
			MPI_Bcast (&current_duration, 1, MPI_INT, 0, local_comm);
			MPI_Bcast (&current_spring, 1, MPI_FLOAT, 0, local_comm);
			


			//md = i % 6 < 5; // MD:VMC = 5:1
			md = rand() > RAND_MAX / 6;
			////////////////////////////////////////////////
			// Effective start of sampling loop           //
			////////////////////////////////////////////////
			MPI_Bcast (&md, 1, MPI_INT, 0, local_comm);

			//Take a sample
			lammps_start_time = get_time();
			if (!md) {
				dx = ((double) rand() / RAND_MAX - 0.5) * 0.01 + 1;
				MPI_Bcast (&dx, 1, MPI_DOUBLE, 0, local_comm);
				sprintf (line, "change_box all x scale %f y scale %f z scale %f", dx, dx, dx);
				lmp->input->one (line);
				lmp->input->one ("change_box all remap");
				U = *((double *) lammps_extract_variable(lmp, (char*)"U", (char*)"all"));
				V = *((double *) lammps_extract_variable(lmp, (char*)"V", (char*)"all"));
			}
			sprintf (line, "run %d", md ? current_duration : 0);
			lmp->input->one (line);
			lammps_split = get_time() - lammps_start_time;
			Q6 = *((double *) lammps_extract_compute(lmp,(char*)"Q6", 0, 0));

			// Compute acceptance
			if (local_rank == 0) {
				d_Q = Q6 - targets[window_index];
				bias_potential_new = d_Q * d_Q;
				log_boltz_factor = (-0.5 * spring_init / p->temperature) * \
								   (bias_potential_new-bias_potential_old);
				accept = log((double) rand() / RAND_MAX) < log_boltz_factor;
			}
			MPI_Bcast (&accept, 1, MPI_INT, 0, local_comm);

			// Write things down
			if (local_rank == 0){
				if (md) {
					fprintf (local_log, "%-9d%-15f%-15f%-4d%-4d%-15lld%-15lld%-8d%-15f\n", \
							i, Q6_old, Q6, accept, !md, lammps_split, \
							get_time() - loop_start_time, current_duration, \
							current_spring);
				} else {
					fprintf (local_log, "%-9d%-15f%-15f%-4d%-4d%-15f%-15f%-15f%15f\n", \
							i, Q6_old, Q6, accept, !md, V_old, V, U_old, U);
				}
			}

			last_step_end_time = get_time();//TODO these got rearranged

			// Unsample if rejected, store state if accepted
			if (!accept) {// Umbrella reject
				if (md) { // Return to last accepted state
					lammps_scatter_atoms(lmp,(char*)"x",1,3,positions_buffer);
					if (local_rank == 0) seed = rand() % 1000 + 1; //TODO why 1000?
					MPI_Bcast (&seed, 1, MPI_INT, 0, local_comm);
					sprintf(line, "velocity all create %f %d", p->temperature, seed);
					lmp->input->one (line);
				} else { // Undo box modulation
					rec = 1 / dx;
					sprintf (line, "change_box all " \
							"x scale %f y scale %f z scale %f", rec, rec, rec);
					lmp->input->one (line);
					lmp->input->one ("change_box all remap");
				}
			} else { //Umbrella accept
				Q6_old = Q6;
				d_Q_old = Q6_old - targets[window_index];
				bias_potential_old = d_Q_old * d_Q_old;
				if (md) {
					lammps_gather_atoms(lmp,(char*)"x",1,3,positions_buffer);
					++accept_count;
				} else {
					++vmc_accept_count;
					U_old = U;
					V_old = V;
				}
			}
			step_time = get_time() - last_step_end_time;
			pthread_mutex_unlock (&mpi_mutex);
			////////////////////////////////////////////////
			// End of sampling loop                       //
			////////////////////////////////////////////////
		}
		delete [] positions_buffer;
		delete lmp;
		if (local_rank == 0)
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
