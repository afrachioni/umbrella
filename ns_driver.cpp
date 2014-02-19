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

// Command line options parser
#include "cl_parser.h"

// Script parser
#include "parser.h"


// Sampling logger
#include "logger.h"


// MPI window splitter
#include "global.h"


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

#if __cplusplus >= 201103L
#define CPP11 1
#include<random>
#else
#warning "No detected support for C++11.  Using std::rand()."
#endif

// print node ID on Kraken
#ifdef KRAKEN
#include <pmi.h>
#endif

#define printmsg(...) if (global->global_rank == 0) fprintf(stdout, __VA_ARGS__);
#define debugmsg(...) if (DEBUG && global->global_rank == 0) fprintf(stdout, __VA_ARGS__);

using namespace LAMMPS_NS;

int main(int narg, char **arg)
{
		MPI_Init(&narg,&arg);
		int me;
		MPI_Comm_rank (MPI_COMM_WORLD, &me);
		if (me == 0) {
			fprintf (stdout, "                       "
					"-------------------------------------------------\n");
			fprintf (stdout, "                       "
					"|  LAMMPS umbrella sampler, version %s  |\n", VERSION);
			fprintf (stdout, "                       "
					"|           Anthony Frachioni,  1-2013          |\n");
			fprintf (stdout, "                       "
					"|            afrachi1@binghamton.edu            |\n");
			fprintf (stdout, "                       "
					"|       Compiled on "__DATE__", "__TIME__"       |\n");
			fprintf (stdout, "                       "
					"-------------------------------------------------\n\n");
		}
		// Parse command line arguments
		CLParser::verbose = me == 0;
		CLParser *p = new CLParser (narg, arg);
		if (p->parse_error) {
			fprintf (stdout, "Parse errors present, exiting...\n");
			MPI_Abort (MPI_COMM_WORLD, 1);//TODO
		}

		// Split MPI_COMM_WORLD into windows
		Global *global = new Global (MPI_COMM_WORLD, p->windows);
		global->split();

		// Setup LAMMPS instance with initial conditions and settings
		debugmsg ("Creating LAMMPSes...\n");
		char *args[] = {(char*)"foo", (char*)"-screen", (char*)"none", \
			(char*)"-log", (char*)"none"};
		LAMMPS *lmp = new LAMMPS(5,args,global->local_comm);

		// Parse input script
		Parser *parser = new Parser (p->script, lmp, global);
		debugmsg ("Processing input script...\n");
		parser->parse();


		// Set up per-window logging (log names currently set here)
		char line[100];
		mkdir ("logs", S_IRWXU);
		sprintf (line, "logs/log_%d.txt", global->window_index);
		Logger *logger = new Logger(line, parser->nparams, (int)global->window_index, \
				global->local_rank, parser->param_ptrs, \
				parser->nsteps, parser->steps);
		logger->init();

		if (p->log_lammps) {
			sprintf (line, "log logs/log_%d.lammps", global->window_index);
			lmp->input->one(line);
		}


		// Execute global window init
		parser->execute_init();

		int natoms = static_cast<int> (lmp->atom->natoms); // Just for fun
		debugmsg ("Number of atoms on root window: %d\n", natoms);

		// Execute per-step initialization blocks
		for (int i = 0; i < parser->nsteps; ++i)
			if ((parser->steps)[i]->probability)
				(parser->steps)[i]->execute_init();

		//Umbrella definitions
		double log_boltzmann;
		int accept;
		UmbrellaStep *chosen_step;
		int steptype;
		float step_rand, accept_rand;

#if CPP11
		std::default_random_engine generator;
		std::uniform_real_distribution<float> distribution (0, 1);
#endif

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

		if (global->global_rank == 0) {
			struct parameter_pointers management_data;
			management_data.window_index = global->window_index;
			management_data.num_active_windows = global->num_windows; 
			management_data.roots_comm = global->roots_comm;
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
		if (global->local_rank == 0) {
			pthread_mutex_lock (&mpi_mutex);
			MPI_Irecv (&msg, 1, message_type, 0, TAG, global->roots_comm, &req);
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
			if (global->local_rank == 0) {
				MPI_Test (&req, &recv_complete, MPI_STATUS_IGNORE);
				if (recv_complete) {
					fprintf (stderr, "Window %d has recieved a message!\n", global->window_index);
					fprintf (stderr, "\tspring_msg: %d\n", msg.spring_message);
					fprintf (stderr, "\tduration_msg: %d\n", msg.duration_message);
					fprintf (stderr, "\tnew_duration: %d\n", msg.new_duration);
					fprintf (stderr, "\tnew_spring: %f\n", msg.new_spring);
					if (msg.spring_message)
						;//current_spring = msg.new_spring;
					if (msg.duration_message)
						;//current_duration = msg.new_duration;
					recv_complete = 0;
					MPI_Irecv (&msg, 1, message_type, 0, TAG, global->roots_comm, &req);
				}
			}
			//bcast to window
			//MPI_Bcast (&current_duration, 1, MPI_INT, 0, global->local_comm);
			//MPI_Bcast (&current_spring, 1, MPI_FLOAT, 0, global->local_comm);
			//
			//--------------------------------------------------------

			////////////////////////////////////////////////
			// Effective start of sampling loop           //
			////////////////////////////////////////////////


			// Increment LAMMPS MC step counter
			sprintf (line, "variable lu_step equal %d", i);
			lmp->input->one (line);


			// Execute any periodic tasks
			for (int j = 0; j < parser->tasks.size(); ++j)
				(parser->tasks)[j]->execute_task(i);

			// Choose a step to execute
			step_rand = (float) rand() / RAND_MAX;
			MPI_Bcast (&step_rand, 1, MPI_INT, 0, global->local_comm);
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
#if CPP11
			accept_rand = distribution(generator);
#else
			accept_rand = (float) rand() / RAND_MAX;
#endif
			accept = log (accept_rand) < log_boltzmann;
			MPI_Bcast (&accept, 1, MPI_INT, 0, global->local_comm);

			// Act accordingly
			if (accept || UmbrellaStep::force_accept) {
				chosen_step->execute_accept();
				for (int j = 0; j < parser->nparams; ++j)
					(parser->param_ptrs)[j]->notify_accepted();
			} else
				chosen_step->execute_reject();

			// Write things down
			logger->step_taken (i, steptype, accept || UmbrellaStep::force_accept || i == 1);//Clean all this up

			if (UmbrellaStep::force_accept)
				UmbrellaStep::force_accept = 0;
			////////////////////////////////////////////////
			// End of sampling loop                       //
			////////////////////////////////////////////////
			pthread_mutex_unlock (&mpi_mutex);
		}
		delete lmp;
		MPI_Finalize();
}
