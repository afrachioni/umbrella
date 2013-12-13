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

#define VERSION "13.12.11.0"
#define DEBUG 1
#define MAX_FNAME_LENGTH 500
#define DUMP_EVERY_STEPS 100

//Command line options parser
#include "cl_parser.h"

// Script parser
#include "parser.h"


// Sampling logger
#include "logger.h"


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
				"|       Compiled on "__DATE__", "__TIME__"       |\n");
		printmsg ("                       "
				"-------------------------------------------------\n");

		//get node ID if on Kraken
		int nid  = -1;
#ifdef KRAKEN
		PMI_CNOS_Get_nid(me, &nid);
#endif

		//parse command line arguments
		CLParser *p = new CLParser (narg, arg); // XXX No CL options in script version
		p->verbose = me == 0;//turn off for DEBUG?
		if (p->parse_error) {
			printmsg ("Parse errors present, exiting...\n");
			MPI_Finalize ();
		}

		//read init file once to determine number of windows
		//FIXME only root process should do this
		int num_windows = 2; // DEBUG


		//---------------------------------------------------------------------
		// MPI window stuff starts here
		//
		//

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
		MPI_Group world_group, roots_group;
		MPI_Comm_group (MPI_COMM_WORLD, &world_group);
		MPI_Group_incl (world_group, num_active_windows, local_roots, \
				&roots_group);
		MPI_Comm roots_comm;
		MPI_Comm_create (MPI_COMM_WORLD, roots_group, &roots_comm);

		//
		//
		// End MPI window stuff
		//---------------------------------------------------------------------




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

		Parser *parser = new Parser ("in.txt", lmp, local_comm, local_rank, \
				roots_comm, num_windows);
		parser->parse();
#if DEBUG	
		// This should happen at runtime, the user might care
		sprintf (line, "log logs/log_%d.lammps", window_index);
		lmp->input->one(line);
#endif


		// Execute global window init
		parser->execute_init();

		int natoms = static_cast<int> (lmp->atom->natoms); // Just for fun
		debugmsg ("Number of atoms: %d\n", natoms);

		// Execute per-step initialization blocks
		for (int i = 0; i < parser->nsteps; ++i)
			if ((parser->steps)[i]->probability)
				(parser->steps)[i]->execute_init();

		//Umbrella definitions
		double log_boltzmann;
		int accept;
		int accept_count = 0;
		int vmc_accept_count = 0;
		UmbrellaStep *chosen_step;
		int steptype;
		float step_rand, accept_rand;

		sprintf (line, "logs/log_%d.txt", window_index);
		Logger *logger = new Logger(line, parser->nparams, (int)window_index, \
				local_rank, parser->param_ptrs, \
				parser->nsteps, parser->steps);
		logger->init();


		// -----------------------------------------------------------
		//  All this should get moved to a separate class
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
		//
		//------------------------------------------------------------



		printmsg ("Samples away!\n\n");
		//Q6_old is most recently accepted Q6
		for (int i = 0; i < p->count; i++) {
			//--------------------------------------------------------
			//
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
			//
			//--------------------------------------------------------
			
			////////////////////////////////////////////////
			// Effective start of sampling loop           //
			////////////////////////////////////////////////


			// Choose a step to execute
			step_rand = (float) rand() / RAND_MAX;
			MPI_Bcast (&step_rand, 1, MPI_INT, 0, local_comm);
			for (steptype = 0; steptype < parser->nsteps; ++steptype) {
				chosen_step = (parser->steps)[steptype];
				if (chosen_step->rand_min <= step_rand && step_rand < chosen_step->rand_max)
					break;
			}

			// Execute step
			chosen_step->execute_step();

			// Add up Boltzmann factors
			log_boltzmann = 0;
			for (int j = 0; j < parser->nparams; ++j)
				log_boltzmann += (parser->param_ptrs)[j]->compute_boltzmann_factor();

			// Compute acceptance
			accept_rand = (float) rand() / RAND_MAX;
			MPI_Bcast (&accept_rand, 1, MPI_INT, 0, local_comm);
			accept = log (accept_rand) < log_boltzmann;

			// Act accordingly
			if (accept)
				chosen_step->execute_accept();
			else
				chosen_step->execute_reject();

			// Write things down
			logger->step_taken (i, steptype, accept);


			////////////////////////////////////////////////
			// End of sampling loop                       //
			////////////////////////////////////////////////
			pthread_mutex_unlock (&mpi_mutex);
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
