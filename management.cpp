#include<string.h>
#include<unistd.h>
#include<pthread.h>

#include<sys/types.h>
#include<sys/stat.h>
#include<fcntl.h>

#define TAG 487
struct parameter_pointers {
	int window_index;
	int num_active_windows;

	MPI_Comm roots_comm;

	pthread_mutex_t *mpi_mutex_ptr;

	MPI_Datatype message_type;
};
struct management_message {
	int duration_message;
	int spring_message;

	int new_duration;
	float new_spring;
};

void *manage (void *params) {

	struct parameter_pointers my_pointers;
	my_pointers = *((struct parameter_pointers *) params);
	fprintf (stdout, "Hello from management thread on root!\n");
	fprintf (stdout, "Syntax: set_spring <window_index> <new_spring>\n");
	fprintf (stdout, "        set_duration <window_index> <new_duration>\n");

	struct management_message msg;

	MPI_Errhandler_set (MPI_COMM_WORLD, MPI_ERRORS_RETURN);
	int tha_rank;
	MPI_Comm_rank (my_pointers.roots_comm, &tha_rank);
	int roots_size;
	MPI_Comm_size (my_pointers.roots_comm, &roots_size);
	fprintf(stdout, "There are %d processes in roots_comm!\n", roots_size);

	FILE *pipe = fopen ("feed", "r");//stdin;

	char line[500];//TODO somehow get MAX_FNAME_LENGTH from ns_driver.cpp
	char first[500];
	char second[500];
	char third[500];

	int window;
	float new_spring;
	int new_duration;

	while(1) {
		msg.duration_message = 0;
		msg.spring_message = 0;
		fgets (line, 500, pipe);
		while (feof(pipe)) {
			usleep(100000);
			clearerr (pipe);
			fgets (line, 500, pipe);
		}
		sscanf (line, "%s %s %s", &first, &second, &third);
		if (strcmp(line, "\n") == 0)
			continue;
		if (strcmp(first, "set_spring") == 0) {
			fprintf (stdout, "set_spring!\n");
			sscanf (second, "%d", &window);
			sscanf (third, "%f", &msg.new_spring);
			msg.spring_message = 1;
		} else if (strcmp(first, "set_duration") == 0) {
			fprintf (stdout, "set_duration!\n");
			sscanf (second, "%d", &window);
			sscanf (third, "%d", &msg.new_duration);
			msg.duration_message = 1;
		} else {
			fprintf (stdout, "Command not recognized: %s\n", first);
			continue;
		}
		if (window < 0 || window > my_pointers.num_active_windows - 1) {
			fprintf (stdout, "Specified window index is outside acceptable bounds.\n");
			continue;
		}
		MPI_Request req;
		pthread_mutex_lock (my_pointers.mpi_mutex_ptr);
		MPI_Isend (&msg, 1, my_pointers.message_type, window, TAG, my_pointers.roots_comm, &req);
		pthread_mutex_unlock (my_pointers.mpi_mutex_ptr);
		line[0] = '\0';
		first[0] = '\0';
		second[0] = '\0';
		third[0] = '\0';
	}
};
