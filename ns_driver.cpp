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
#include "cl_parser.h"

// Script parser
#include "parser.h"


// LAMMPS include files
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <mpi.h>
#include <lammps.h>
#include <input.h>
#include <atom.h>
#include <library.h>
#include "management.cpp"

// system include files
#include <math.h>
#include <pthread.h>
#include <sys/time.h>
#include <sys/stat.h> //mkdir needs this
#include <sys/types.h> //and this too

// print node ID on Kraken
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
		CLParser *p;
		p = new CLParser (narg, arg);

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
		//p = new CLParser (narg, arg); XXX No CL options in script version
		//if (p->parse_error) {
			//printmsg ("Parse errors present, exiting...\n");
			//MPI_Finalize ();
		//}

		//read init file once to determine number of windows
		//FIXME only root process should do this
		int num_windows = 2; // DEBUG
		/* XXX Keep this code for files on which <<>> is operated
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
		*/
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
		} else {
			printmsg ("There are %d available processors and %d defined "
					"windows; each window will use %d processors.\n", \
					nprocs, num_windows, nprocs / num_windows);
		}

		int num_active_windows;
		if (nprocs < num_windows)
			num_active_windows = nprocs;
		else
			num_active_windows = num_windows;

		// Split COMM_WORLD into communicators for each window
		debugmsg ("At line %d\n", __LINE__);
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

		// Create roots_comm, a communicator for window roots
		debugmsg ("At line %d\n", __LINE__);
		MPI_Group world_group, roots_group;
		MPI_Comm_group (MPI_COMM_WORLD, &world_group);
		MPI_Group_incl (world_group, num_active_windows, local_roots, \
				&roots_group);
		MPI_Comm roots_comm;
		MPI_Comm_create (MPI_COMM_WORLD, roots_group, &roots_comm);




		// Setup LAMMPS instance with initial conditions and settings
		char line[100];
		debugmsg ("Creating LAMMPSes...\n");
		sprintf (line, "logs/log_%d.screen", window_index);
		//char *args[] = {(char*)"foo", (char*)"-screen", \
		line, (char*)"-log", (char*)"none"};
		char *args[] = {(char*)"foo", (char*)"-screen", (char*)"none", \
			(char*)"-log", (char*)"none"};
		LAMMPS *lmp = new LAMMPS(5,args,local_comm);
		//LAMMPS *lmp = new LAMMPS(0,NULL,local_comm);//for all stdout

		Parser *parser = new Parser ("in.txt", lmp);
		parser->parse();
#if DEBUG	
		// This should happen at runtime, the user might care
		sprintf (line, "log logs/log_%d.lammps", window_index);
		lmp->input->one(line);
#endif
		parser->execute_init();

		int natoms = static_cast<int> (lmp->atom->natoms);
		debugmsg ("Number of atoms: %d\n", natoms);


		// Store bounds in step objects?
		for (int i = 0; i < parser->nsteps; ++i)
			fprintf (stderr, "P[%d]: %f\n", (parser->steps)[i]->probability);
		//Umbrella definitions
		int accept = 1;
		int seed;
		int accept_count = 0;
		int vmc_accept_count = 0;
		int64_t last_step_end_time, step_time, lammps_start_time, lammps_split;
		last_step_end_time = get_time();
		//double Q6_old = Q6;
		double d_Q, d_Q_old, bias_potential_new, bias_potential_old;

		const int64_t loop_start_time = get_time();


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

			const int update_count = 4;
			int update_lengths[update_count] = {1, 1, 1, 1};
			MPI_Aint update_offsets[update_count] = {0, sizeof(int), 2*sizeof(int), 2*sizeof(int) + sizeof(double)};
			MPI_Datatype update_types[update_count] = {MPI_INT, MPI_INT, MPI_DOUBLE, MPI_DOUBLE};
			MPI_Datatype update_message_type;
			MPI_Type_struct (update_count, update_lengths, update_offsets, update_types, &update_message_type);
			MPI_Type_commit (&update_message_type);

		//mkfifo ("logs/feed_pipe", S_IWUSR | S_IRUSR);
		//int fd = open ("logs/feed_pipe", O_WRONLY | O_NONBLOCK);



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
						;//current_spring = msg.new_spring;
					if (msg.duration_message)
						;//current_duration = msg.new_duration;
					recv_complete = 0;
					MPI_Irecv (&msg, 1, message_type, 0, TAG, roots_comm, &req);
				}
			}
			//bcast to window
			//MPI_Bcast (&current_duration, 1, MPI_INT, 0, local_comm);
			//MPI_Bcast (&current_spring, 1, MPI_FLOAT, 0, local_comm);
			
			////////////////////////////////////////////////
			// Effective start of sampling loop           //
			////////////////////////////////////////////////


			step_time = get_time() - last_step_end_time;
			pthread_mutex_unlock (&mpi_mutex);
			////////////////////////////////////////////////
			// End of sampling loop                       //
			////////////////////////////////////////////////
		}
		delete lmp;
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
