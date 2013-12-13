#ifndef WINDOW_INFO_H
#define WINDOW_INFO_H

#include <mpi.h>

class Global {
	public:
		Global (MPI_Comm world, int num_windows);
		int nprocs;
		int global_rank;
		int local_rank;
		int num_windows;
		int window_index;

		MPI_Comm world;
		MPI_Comm local_comm;
		MPI_Comm roots_comm;

		void split();
		void abort(char *message);
		void warn(char *message);
};
#endif
