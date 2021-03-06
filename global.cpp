#include "global.h"

#include <stdio.h>
#include <stdlib.h>
#include <cstring>

#define printmsg(...) if (global_rank == 0) fprintf(stdout, __VA_ARGS__);

// static
Global *Global::instance = NULL;

// static
Global* Global::get_instance() {
	if (instance == NULL)
		fprintf (stderr,
		"############################################################"
		"\nget_instance called on Global without first creating one."
		"\nThis is an internal error which should never happen.\n"
		"############################################################\n");
	return instance;
}

// static
void Global::init (MPI_Comm world, int num_windows) {
	instance = new Global (world, num_windows);
	instance->split();
}

Global::Global (MPI_Comm world, int num_windows) {
	this->world = world;
	this->num_windows = num_windows;
	MPI_Win_create (&abort_called, sizeof(int), sizeof(int), \
				MPI_INFO_NULL, world, &window);
	abort_called = 0;
}

void Global::finalize () {
	MPI_Win_free (&window);
	MPI_Finalize();
}

int Global::get_global_rank() {
	return global_rank;
}

int Global::get_local_rank() {
	return local_rank;
}

int Global::get_num_windows() {
	return num_windows;
}

int Global::get_window_index() {
	return window_index;
}

void Global::split() {
	MPI_Comm_rank(world, &global_rank);
	MPI_Comm_size(world, &nprocs);
	char line[200];
	if (nprocs < num_windows) {
		sprintf (line, "More windows have been defined (%d) than "
				"exist available processors (%d).\n\t\tThe first %d windows "
				"will execute.", num_windows, nprocs, nprocs);
		warn (line);
	} else if (nprocs % num_windows != 0) {
		sprintf (line, "The number of available processors (%d) is "
				"not evenly divisible by the number of defined windows "
				"(%d).\n\t\tWindows 0 to %d will use %d "
				"processors, windows %d to %d will use %d "
				"processors.", nprocs, num_windows, \
				nprocs % num_windows - 1, nprocs / num_windows + 1, \
				nprocs % num_windows, num_windows - 1, \
				nprocs / num_windows);
		warn (line);
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
	MPI_Comm_split(world, global_rank % num_active_windows, global_rank, &local_comm);
	MPI_Comm_rank (local_comm, &local_rank);
	window_index = global_rank % num_active_windows;
	int local_roots_send[num_active_windows];
	int local_roots[num_active_windows];
	for (int i = 0; i < num_active_windows; ++i)
		local_roots_send[i] = 0;
	if (local_rank == 0)
		local_roots_send[window_index] = global_rank;
	MPI_Allreduce (local_roots_send, local_roots, num_active_windows, \
			MPI_INT, MPI_SUM, world);

	// Create roots_comm, a communicator for window roots
	MPI_Group world_group, roots_group;
	MPI_Comm_group (world, &world_group);
	MPI_Group_incl (world_group, num_active_windows, local_roots, \
			&roots_group);
	MPI_Comm_create (world, roots_group, &roots_comm);
}

// Halt entire job gracefully.  Must be called by all processes.
void Global::stop(char *message) {
	if (global_rank == 0)
		fprintf (stdout, "\033[31m%s\033[0m\n", message);
	finalize();
	exit(1);
}

// Fancy wrapper to MPI_Abort.  Can be called by any process, and will only
// print message from first caller, which then solely calls MPI_Abort.
void Global::abort(char *message) {
	if (global_rank) {
		MPI_Win_lock (MPI_LOCK_EXCLUSIVE, 0, 0, window);
		MPI_Get (&abort_called, 1, MPI_INT, 0, 0, 1, MPI_INT, window);
		if (!abort_called) {
			int one = 1;
			MPI_Put (&one, 1, MPI_INT, 0, 0, 1, MPI_INT, window);
		}
		MPI_Win_unlock (0, window);
	}
	if (abort_called) MPI_Barrier (world); // all going to die.
	else {
		fprintf (stderr, "\033[31m\n");
		// TODO N be safe
		char line[500];
		for (unsigned i = 0; i < strlen(message)+13; ++i)
			line[i] = '=';
		line[strlen(message)+13] = '\n';
		line[strlen(message)+14] = '\0';
		fputs (line, stderr);
		fprintf (stderr, "|| ERROR: %s ||\n", message);
		fputs (line, stderr);
		fprintf (stderr, "\033[0m\n");
		fprintf (stderr, "\nKilling %d processes...\n\n", nprocs);
		MPI_Abort (world, 1);
	}
}

void Global::warn(char *message) {
	if (global_rank == 0)
		fprintf (stdout, "\033[33mWARNING: %s\033[0m\n", message);
}

void Global::debug (char *message) {
	if (global_rank == 0)
		fprintf (stderr, "\033[32m%s\033[0m\n", message);
}
