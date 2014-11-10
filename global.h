#ifndef GLOBAL_H
#define GLOBAL_H

#include <mpi.h>

// This is a singleton, mostly to protect me from my own future stupidity.
// Fields cannot be accessed without a copy of the instance, the interface to 
// which is NULL until the init method is called.

class Global {
	public:
		static Global *get_instance();
		static void init (MPI_Comm world, int num_windows);

		// TODO maybe make these static?
		int get_global_rank();
		int get_local_rank();
		int get_num_windows();
		int get_window_index();

		MPI_Comm world;
		MPI_Comm local_comm;
		MPI_Comm roots_comm;

		void split();
		void stop(char *message);
		void abort(char *message);
		void warn(char *message);
		void debug (char *message);
		void finalize ();
	private:
		int nprocs;
		int global_rank;
		int local_rank;
		int num_windows;
		int window_index;
		static Global *instance;
		int abort_called;
		MPI_Win window;
		Global (MPI_Comm world, int num_windows);
		Global (Global const&); // Make no copies! (Not implemented.)
		void operator = (Global const&); // Or assignments.
};
#endif
