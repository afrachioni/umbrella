#include "global.h"

#include <stdio.h>

#define printmsg(...) if (global_rank == 0) fprintf(stdout, __VA_ARGS__);

Global::Global (MPI_Comm world, int num_windows) {
	this->world = world;
	this->num_windows = num_windows;
}

void Global::split() {
		MPI_Comm_rank(MPI_COMM_WORLD,&global_rank);
		MPI_Comm_size(MPI_COMM_WORLD,&nprocs);
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
		MPI_Comm_split(MPI_COMM_WORLD, global_rank % num_active_windows, global_rank, &local_comm);
		MPI_Comm_rank (local_comm, &local_rank);
		window_index = global_rank % num_active_windows;
		int local_roots_send[num_active_windows];
		int local_roots[num_active_windows];
		for (int i = 0; i < num_active_windows; ++i)
			local_roots_send[i] = 0;
		if (local_rank == 0)
			local_roots_send[window_index] = global_rank;
			MPI_Allreduce (local_roots_send, local_roots, num_active_windows, \
					MPI_INT, MPI_SUM, MPI_COMM_WORLD);

		// Create roots_comm, a communicator for window roots
		MPI_Group world_group, roots_group;
		MPI_Comm_group (world, &world_group);
		MPI_Group_incl (world_group, num_active_windows, local_roots, \
				&roots_group);
		MPI_Comm_create (MPI_COMM_WORLD, roots_group, &roots_comm);
}

void Global::abort(char *message) {
	if (global_rank == 0) {
		fprintf (stderr, "\033[31m\n");
		fprintf (stderr, \
		"=====================================================================\n");
		fprintf (stderr, "ERROR: %s\n", message);
		fprintf (stderr, \
		"=====================================================================\n");
		fprintf (stderr, "\033[0m\n");
		fprintf (stderr, "\nKilling %d processes...\n\n", nprocs);
	}
	MPI_Barrier (MPI_COMM_WORLD);
	MPI_Abort (MPI_COMM_WORLD, 1);
}

void Global::warn(char *message) {
	if (global_rank == 0)
		fprintf (stdout, "\033[33mWARNING: %s\033[0m\n", message);
}
