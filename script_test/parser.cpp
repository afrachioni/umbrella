#include <stdio.h>
#include <string.h>
#include <vector>
#include <map>
#include <string>

#include "umbrella_step.h"
#include "umbrella_parameter.h"
#include "parser.h"

#define MAX_LINE_LENGTH 1000


int main (int nargs, char **args) {
	MPI_Init(&nargs,&args);
	Parser *p = new Parser("in.txt", NULL);
	p->parse();
	int me;
	MPI_Comm_rank(MPI_COMM_WORLD, &me);
	if (me == 0)
		p->print();
	delete p;
	MPI_Finalize();
}

Parser::Parser(const char *fname, LAMMPS_NS::LAMMPS *) {
	strcpy (this->fname, fname);
	this->lmp = lmp;
};

Parser::~Parser() {
	// XXX can't figure out why this throws a segmentation fault./
	//std::vector<std::string> *s;
	//for (std::map<std::string, UmbrellaStep>::iterator it = steps_map.begin(); it != steps_map.end(); ++it) {
		//fprintf (stderr, "Deleting UmbrellaStep: %s\n", it->first.c_str());
		//delete &(it->second);
	//}
		//delete &steps_map;
};

void Parser::parse() {
	FILE *in_p = fopen (fname, "r");
	int n;
	int length;
	char line_buf[MAX_LINE_LENGTH];
	char *line;
	char first_token[20];
	char second_token[20];
	char third_token[20];
	char fourth_token[20];
	char fifth_token[20];
	UmbrellaParameter *p;

	int me;
	MPI_Comm_rank(MPI_COMM_WORLD, &me);
	char *file_data;
	int max_line_length = 0;
	int num_lines;

	// Root reads file, everybody parses (because passing structs is annoying)
	std::vector<std::string> line_ptrs;
	if (me == 0) {
		// Read file into vector<string>, keeping track of longest line
		int line_length;
		while (1) {
			fgets (line_buf, MAX_LINE_LENGTH, in_p);
			if (feof (in_p)) break;
			line_length = strlen (line_buf);
			line_buf[line_length - 1] = '\0';
			std::string m(line_buf);
			line_ptrs.push_back (m);
			if (line_length > max_line_length)
				max_line_length = line_length;
		}
		fclose (in_p);
		num_lines = line_ptrs.size();
	}
	// Broadcast to all my friends
	MPI_Bcast ( &max_line_length, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast ( &num_lines, 1, MPI_INT, 0, MPI_COMM_WORLD);
	// Copy to contiguous buffer for shipment over the network
	file_data = new char [max_line_length * num_lines];
	if (me == 0)
		for (int i = 0; i < line_ptrs.size(); ++i)
			strcpy (& file_data [i * max_line_length], line_ptrs[i].c_str());
	MPI_Bcast ( file_data, max_line_length * num_lines, MPI_CHAR, 0, MPI_COMM_WORLD);

	std::vector<std::string>* current_block = NULL;
	// no exceptions yet for too few tokens or things not specified
	// Loop over lines in file buffer, populate appropriate structures
	for (int i = 0; i < num_lines; ++i) {
			line = file_data + i * max_line_length;
			length = strlen (line);
			n = sscanf (line, "%s %s %s %s", first_token, second_token, third_token, fourth_token);
		if (strcmp (first_token, "#AF") == 0 && n > 0) {
			//fprintf (stderr, "Special directive: %s\n", line);
			if (n == 1) {
				fprintf (stderr, "Parse error: empty directive at line %d.\n", i);
				break;
			} else if (strcmp (second_token, "steptype") == 0) {
				UmbrellaStep *s = new UmbrellaStep (lmp); // TODO when does this die?
				steps_map[third_token] = *s;

			} else if (strcmp (second_token, "parameter") == 0) {
				p = new UmbrellaParameter (third_token, fourth_token, fifth_token, lmp);
				params.push_back (*p);

			} else if (strcmp (second_token, "takestep") == 0) {
				if (steps_map.find(third_token) == steps_map.end()) {
					fprintf (stderr, "Parse error: takestep before %s defined\n", third_token);
					break;
				}
				current_block = steps_map[third_token].get_take_step_block();
			} else if (strcmp (second_token, "ifaccept") == 0) {
				if (steps_map.find(third_token) == steps_map.end()) {
					fprintf (stderr, "Parse error: ifaccept before %s defined\n", third_token);
					break;
				}
				current_block = steps_map[third_token].get_if_accept_block();
			} else if (strcmp (second_token, "ifreject") == 0) {
				if (steps_map.find(third_token) == steps_map.end()) {
					fprintf (stderr, "Parse error: ifreject before %s defined\n", third_token);
					break;
				}
				current_block = steps_map[third_token].get_if_reject_block();
			} else if (strcmp (second_token, "stepinit") == 0) {
				if (steps_map.find(third_token) == steps_map.end()) {
					fprintf (stderr, "Parse error: stepinit before %s defined\n", third_token);
					break;
				}
				current_block = steps_map[third_token].get_step_init_block();
			} else {
				fprintf (stderr, "Directive not recognized: %s\n", line);
			}
		}
		if (line[0] == '#' || line[0] == '\0') continue; //perhaps pass to LAMMPS so they show up on the logs
		if (current_block != NULL) {
			current_block->push_back (line);
		}
	}
	delete [] file_data;
}

void Parser::print() {
	fprintf (stdout, "\n\n\nDefined steps:\n");
	std::vector<std::string> *v;
	for (std::map<std::string, UmbrellaStep>::iterator it = steps_map.begin(); it != steps_map.end(); ++it) {
		fprintf ( stdout, "\tStep type: %s\n", it->first.c_str());
		v = it->second.get_take_step_block();
		fprintf (stdout, "\t\tTake step:\n");
		for (int j = 0; j < v->size(); ++j)
			fprintf (stdout, "\t\t\t%s\n", (*v)[j].c_str());//v->at(j));

		v = it->second.get_if_reject_block();
		fprintf (stdout, "\t\tIf reject:\n");
		for (int j = 0; j < v->size(); ++j)
			fprintf (stdout, "\t\t\t%s\n", v->at(j).c_str());

		v = it->second.get_if_accept_block();
		fprintf (stdout, "\t\tIf accept:\n");
		for (int j = 0; j < v->size(); ++j)
			fprintf (stdout, "\t\t\t%s\n", v->at(j).c_str());

		v = it->second.get_step_init_block();
		fprintf (stdout, "\t\tStep init:\n");
		for (int j = 0; j < v->size(); ++j)
			fprintf (stdout, "\t\t\t%s\n", v->at(j).c_str());
	}
	fprintf (stdout, "Defined parameters:\n");
	for (int j = 0; j < params.size(); ++j)
		fprintf (stdout, "\t%s\n", j, params[j].param_vname);
}
