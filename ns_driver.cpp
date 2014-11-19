/* ----------------------------------------------------------------------
                     LAMMPS umbrella sampling driver
   Links against a modified LAMMPS to execute umbrella sampling
   Uses custom "compute boop" to calculate order parameter in situ
   Built to examine nucleation in tin

   Copyright 2013 Anthony Frachioni - All Rights Reserved
                          afrachioni@gmail.com
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

// Command line options parser
#include "cl_parser.h"

// Script parser
#include "parser.h"


// Sampling logger
#include "logger.h"


// MPI window splitter
#include "global.h"

// For final report
#include "barostat_step.h"

// Histogram
#include "histogram.h"

// LAMMPS include files
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <mpi.h>
#include <lammps.h>
#include <input.h>
#include <atom.h>
#include <library.h>
#include <version.h>
//#include "management.cpp"

// system include files
#include <math.h>
#include <pthread.h>
#include <sys/time.h>
#include <sys/stat.h> //mkdir needs this
#include <sys/types.h> //and this too

#define printmsg(...) if (global->get_global_rank() == 0) fprintf(stdout, __VA_ARGS__);
#define debugmsg(...) if (DEBUG && global->get_global_rank() == 0) fprintf(stdout, __VA_ARGS__);

using namespace LAMMPS_NS;

extern const char *gitversion;
extern const char *gitmessage;
extern const char *gitbranch;

int main(int narg, char **arg)
{
		MPI_Init(&narg,&arg);
		int me;
		MPI_Comm_rank (MPI_COMM_WORLD, &me);
		if (me == 0) {
			fprintf (stdout, "                       "
					"-------------------------------------------------\n");
			fprintf (stdout, "                       "
					"|            LAMMPS umbrella sampler            |\n");
			fprintf (stdout, "                       "
					"|           Anthony Frachioni, 1-2013           |\n");
			fprintf (stdout, "                       "
					"|              afrachioni@gmail.com             |\n");
			fprintf (stdout, "                       "
					"|                                               |\n");
			fprintf (stdout, "                       "
					"| LAMMPS version: %-29s |\n", LAMMPS_VERSION);
			fprintf (stdout, "                       "
					"| Latest commit: %-30s |\n", gitmessage);
			fprintf (stdout, "                       "
					"| On branch: %-34s |\n", gitbranch);
			fprintf (stdout, "                       "
					"|    %s   |\n", gitversion);
			fprintf (stdout, "                       "
					"-------------------------------------------------\n\n");
		}
		// Parse command line arguments
		CLParser::verbose = me == 0;
		CLParser *p = new CLParser (narg, arg);
		if (p->parse_error) {
			if (me == 0)
				fprintf (stdout, "Command line errors present, exiting...\n");
			MPI_Finalize ();
			return 0;
		}

		// Split MPI_COMM_WORLD into windows
		Global::init (MPI_COMM_WORLD, p->windows);
		Global *global = Global::get_instance();
		srand(me);

		// Warn about hardcoded integrate accept/reject as global
		//global->warn("Temperature hardcoded to 10K");
		global->warn((char*)"Step named \"integrate\" hardcoded to provide global"
				" accept/reject blocks for now");
		//global->warn("Using Lennard-Jones reduced units!  (Compiled in.)");
		//global->warn("Hardcoded to never bias. (Ever.)  (Really.)");

		// Setup LAMMPS instance with initial conditions and settings
		debugmsg ("Creating LAMMPSes...\n");

		// Original -- wierd numbers in logs
		char *args[] = {(char*)"foo", (char*)"-screen", (char*)"none", \
			(char*)"-log", (char*)"none"};
		LAMMPS *lmp = new LAMMPS(5,args,global->local_comm);

		
		/*
		 * sterilize this just in case
		// Also seems to work
		//char *args[] = {(char*)"foo", (char*)"-echo", (char*)"none", \
			(char*)"-log", (char*)"none"};
		//LAMMPS *lmp = new LAMMPS(5,args,global->local_comm);

		// Seems to work
		//char *args[] = {(char*)"foo", \
			(char*)"-log", (char*)"none"};
		//LAMMPS *lmp = new LAMMPS(3,args,global->local_comm);

		*/
		// Parse input script
		Parser *parser = new Parser (p->script, lmp);
		debugmsg ("Processing input script...\n");
		if (parser->parse()) {
			if (me == 0)
				fprintf (stdout, "\nScript errors present:\n");
			global->stop(parser->error_message());
		}


		// Set up per-window logging (log names currently set here)
		char line[100];
		mkdir ("logs", S_IRWXU);
		sprintf (line, "logs/log_%d.txt", global->get_window_index());
		Logger *logger = new Logger(line, parser->nparams, (int)global->get_window_index(), \
				global->get_local_rank(), parser->param_ptrs, \
				parser->nsteps, parser->steps);
		logger->init();

		if (p->log_lammps) {
			sprintf (line, "log logs/log_%d.lammps", global->get_window_index());
			lmp->input->one(line);
		}

		// Execute global window init
		parser->execute_init();

		int natoms = static_cast<int> (lmp->atom->natoms); // Just for fun
		debugmsg ("Number of atoms on root window: %d\n", natoms);

		// TODO maybe this belongs somewhere else
		sprintf (line, "variable lu_natoms equal %d", natoms);
		lmp->input->one (line);
		lmp->input->one ("variable lu_vol equal vol");

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

		// -----------------------------------------------------------
		//  All this should get moved to a separate class
		//pthread_mutex_t mpi_mutex = PTHREAD_MUTEX_INITIALIZER;
		/*

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
			//if (pthread_create( &manager, NULL, &manage, &management_data))
				//printmsg ("Could not create management thread on root!");
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
		*/
		//
		//------------------------------------------------------------

		int local_count = 0;
		int local_accept_count = 0;


		mkdir ("hist_data", S_IRWXU);
		sprintf (line, "hist_data/window_%d", global->get_window_index());
		mkdir (line, S_IRWXU);
		mkdir ("window_stats", S_IRWXU);
		sprintf (line, "window_stats/window_%d", global->get_window_index());
		mkdir (line, S_IRWXU);

		printmsg ("Samples away!\n\n");
		//Q6_old is most recently accepted Q6
		int64_t start_time = Logger::get_time();
		int64_t step_start_time = Logger::get_time();
		for (int i = 0; i < p->count + 1; ++i) {
			if (i % 1 == 0 && global->get_global_rank() == 0) {
				int64_t now = Logger::get_time();
				int64_t split = now - step_start_time;
				step_start_time = now;
				double rate = local_count?(double)local_accept_count/local_count:-1;
				double barostat_rate = BarostatStep::get_rate();
				BarostatStep::zero_rate();
				fprintf (stdout, "Step: %d\tRate: %f\tBarostat rate: %f\tSplit: %" PRId64 "\n", i, rate, barostat_rate, split);
				local_accept_count = 0;
				local_count = 0;
				//global->debug (line);
			}
			//--------------------------------------------------------
			//
			/* Did this recently due to wierd things
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
			*/
			//bcast to window
			//MPI_Bcast (&current_duration, 1, MPI_INT, 0, global->local_comm);
			//MPI_Bcast (&current_spring, 1, MPI_FLOAT, 0, global->local_comm);
			//
			//--------------------------------------------------------

			//|\|/|\|/|\|/|\|/|\|/|\|/|\|/|\|/|\|/|\|/|\|/|\|
			// Effective start of sampling loop             |
			//|\|/|\|/|\|/|\|/|\|/|\|/|\|/|\|/|\|/|\|/|\|/|\|


			// Increment LAMMPS MC step counter
			sprintf (line, "variable lu_step equal %d", i);
			lmp->input->one (line);


			// Execute any periodic tasks
			for (unsigned j = 0; j < parser->tasks.size(); ++j)
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

			// Skip bias if necessary
			if (i % parser->bias_every || 0)
				continue;

			// Add up Boltzmann factors
			log_boltzmann = 0;
			for (int j = 0; j < parser->nparams; ++j)
				log_boltzmann += (parser->param_ptrs)[j]->compute_boltzmann_factor();

			// Compute acceptance
			accept_rand = (float) rand() / RAND_MAX;
			accept = log (accept_rand) < log_boltzmann;
			MPI_Bcast (&accept, 1, MPI_INT, 0, global->local_comm);

			// Act accordingly
			++local_count;
			if (accept || UmbrellaStep::force_accept) {
				++local_accept_count;
				//global->debug ("\t\t\t\tACCEPT");
				if (UmbrellaStep::force_accept)
					global->debug((char*)"||||||||||||||FORCE||||||||||||||");
				//chosen_step->execute_accept();  // TODO switch this to a global accept iff bias_every > 1

				// XXX Debug, experimental
				parser->steps_map["integrate"]->execute_accept();
				for (int j = 0; j < parser->nparams; ++j) {
					(parser->param_ptrs)[j]->notify_accepted();
				}
			} else {
				//chosen_step->execute_reject();  // TODO switch this to a global reject iff bias_every > 1

				// XXX Debug, experimental
				parser->steps_map["integrate"]->execute_reject();
			}

			// Write things down
			logger->step_taken (i, steptype, accept || UmbrellaStep::force_accept);

			for (std::vector<Histogram *>::iterator it = parser->histograms.begin();
					it != parser->histograms.end(); ++it) {
				(*it)->update();
				(*it)->write(i); //TODO make this better
			}

			if (UmbrellaStep::force_accept)
				UmbrellaStep::force_accept = 0;
		}
			//
			//
			//
			//
			////////////////////////////////////////////////
			// End of sampling loop                       //
			////////////////////////////////////////////////

			//pthread_mutex_unlock (&mpi_mutex);

		delete lmp; // Needs to be alive for extraction!
		delete parser;

		if (global->get_global_rank() == 0) {
			int64_t walltime = Logger::get_time() - start_time;
			fprintf (stdout, "Root window has finished sampling.\n");
			fprintf (stdout, "Walltime / s: %" PRId64 "\n", walltime/1000);
			fprintf (stdout, "Mean sample period / s: %e\n", \
				walltime/((float) p->count*1000));
			fprintf (stdout, "Acceptance rate among volume moves: %f\n",
					BarostatStep::get_rate());
			fprintf (stdout, "Waiting for other windows...\n");
		}
		Global::get_instance()->finalize();
		return 0;
}
