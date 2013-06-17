#include<string.h>

struct parameter_pointers {
	int *duration;
	float *spring;

	int rank;
	int local_rank;
	int window_index;
	int num_active_windows;
	int *local_roots;

	MPI_Comm roots_comm;
	MPI_Comm local_comm;
};
struct management_message {
	int duration_message;
	int spring_message;

	int new_duration;
	int new_spring;
};

void *manage (void *params) {
	const int count = 4;
	int lengths[count] = {1, 1, 1, 1};
	MPI_Aint offsets[count] = {0, sizeof(int), 2*sizeof(int), 3*sizeof(int)};
	MPI_Datatype types[count] = {MPI_INT, MPI_INT, MPI_INT, MPI_INT};
	MPI_Datatype message_type;
	MPI_Type_struct (count, lengths, offsets, types, &message_type);
	MPI_Type_commit (&message_type);



	struct parameter_pointers my_pointers;
	my_pointers = *((struct parameter_pointers *) params);

	struct management_message msg;

	MPI_Errhandler_set (MPI_COMM_WORLD, MPI_ERRORS_RETURN);
	

	while(1) {
		//If global root
		if (my_pointers.rank == 0) {
			fprintf (stdout, "Hello from management thread on root!\n");
			fprintf (stdout, "Syntax: set_spring <window_index> <new_spring>\n");
			fprintf (stdout, "        set_duration <window_index> <new_duration>\n");
			char line[500];//TODO somehow get MAX_FNAME_LENGTH from ns_driver.cpp
			char first[500];
			char second[500];
			char third[500];
			msg.duration_message = 0;
			msg.spring_message = 0;
			int window;
			float new_spring;
			int new_duration;
			fgets (line, 500, stdin);
			sscanf (line, "%s %s %s", &first, &second, &third);
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
			if (window != 0) {
				int err;
				fprintf (stdout, "Process %d right before send\n", my_pointers.rank);
				err = MPI_Send (&msg, 1, message_type, my_pointers.local_roots[window], 0, \
						MPI_COMM_WORLD);
				if (err != MPI_SUCCESS) {
					char err_string[500];
					int err_string_length;
					MPI_Error_string (err, err_string, &err_string_length);
					fprintf (stderr, "MPI ERROR [%d]: %s\n", my_pointers.rank, err_string);
				} else
					fprintf (stdout, "Message sent by process of rank %d!\n", my_pointers.rank);
			}
		} else if (my_pointers.local_rank == 0 && my_pointers.rank != 0) {
			MPI_Status status;
			int err;
			fprintf (stdout, "Process %d right before recv\n", my_pointers.rank);
			MPI_Errhandler_set (MPI_COMM_WORLD, MPI_ERRORS_RETURN);
			err = MPI_Recv (&msg, 1, message_type, 0, 0, MPI_COMM_WORLD, &status);
			if (err != MPI_SUCCESS) {
				char err_string[500];
				int err_string_length;
				MPI_Error_string (err, err_string, &err_string_length);
				fprintf (stderr, "MPI ERROR [%d]: %s\n", my_pointers.rank, err_string);
			} else
				fprintf (stdout, "Message recieved by process of rank %d!\n", my_pointers.rank);
		}
		/*
		fprintf (stdout, "Process %d right before broadcast\n", my_pointers.rank);
		MPI_Bcast (&msg, 1, message_type, 0, my_pointers.local_comm);
		fprintf (stdout, "Process %d right after broadcast\n", my_pointers.rank);
		if (msg.spring_message)
			*(my_pointers.spring) = msg.new_spring;
		if (msg.duration_message)
			*(my_pointers.duration) = msg.new_duration;
			*/
	}
};
//non-root workers block on bcast while their local root blocks on recv
