#include <iostream>

#include <stdio.h>
#include <stdlib.h> //strtod
#include <cstdlib>
#include <string.h>
#include <vector>
#include <map>
#include <string>

#include "parser.h"
#include "quantity.h"
#include "barostat_step.h"
#include "global.h"


#define MAX_LINE_LENGTH 1000
#define MAX_TOKENS 100


/*
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
*/

Parser::Parser(const char *fname, LAMMPS_NS::LAMMPS *lmp) {
	strcpy (this->fname, fname);
	this->lmp = lmp;
	msg[0] = '\0';
	bias_every = 1;
};

Parser::~Parser() {
	std::map<std::string, UmbrellaStep*>::iterator it;
	for (it = steps_map.begin(); it != steps_map.end(); ++it)
		delete (it->second);
};

int Parser::parse() {
	FILE *in_p = fopen (fname, "r");
	if (in_p == NULL) {
		char msg[100];
		sprintf (msg, "Cannot read from file: %s", fname);
		Global::get_instance()->abort(msg);
	}
	char line_buf[MAX_LINE_LENGTH];


	int me;
	MPI_Comm_rank(MPI_COMM_WORLD, &me);
	int max_line_length = 0;
	int num_lines;

	// Root reads file, everybody parses
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
	char *file_data = new char [max_line_length * num_lines]; //TODO delete
	if (me == 0)
		for (unsigned i = 0; i < line_ptrs.size(); ++i)
			strcpy (& file_data [i * max_line_length], line_ptrs[i].c_str());
	MPI_Bcast ( file_data, max_line_length * num_lines,
					MPI_CHAR, 0, MPI_COMM_WORLD);

	std::vector<std::string>* current_block = &init_block;
	// TODO no exceptions yet for too few tokens or things not specified
	// Loop over lines in file buffer, populate appropriate structures
	char *e; // error for string to number conversion
	char *line;
	Quantity *temp = new Quantity((char*)"c_thermo_temp", lmp, true, false);
	for (int i = 0; i < num_lines; ++i) {
		int ln = i + 1;
		line = file_data + i * max_line_length;

		std::string sx = Parser::process_brackets (line);
		line = (char*) sx.c_str();

		strcpy(line_buf, line);

		int n = 0;
		char **tokens = new char*[MAX_TOKENS];//TODO delete
		tokens[n] = strtok(line_buf, " ");
		while (tokens[n] != NULL) {
				tokens[++n] = strtok(NULL, " ");
		}

		char *term = new char ('\0');
		for (int j = 0; j < 100; ++j)
				if (tokens[j] == NULL)
						tokens[j] = term;

		if (n > 0 && strcmp (tokens[0], "#AF") == 0) {
			if (n == 1) {
				sprintf(msg, "Parse error: empty directive at line %d.\n", ln);
				return 1;
			} else if (strcmp (tokens[1], "temperature") == 0) {
				if (params.size() > 0) {
					fprintf (stdout, "Parse warning: non-default temperature "
							"used for some, but not all, parameters");
				}
				temp = new Quantity (tokens[2], lmp, true, false);
				// TODO format all parse errors like this
				if (!temp->is_valid())
					sprintf(msg, "%s\nerror: \"%s\" is not a valid temperature"
							" (line %d).", msg, tokens[2], ln);
			} else if (strcmp (tokens[1], "step_type") == 0) {
				// TODO complain about number of tokens if necessary
				Quantity *d = new Quantity (tokens[3], lmp, true, false);
				if (!d->is_valid() || !d->is_constant())
					sprintf (msg, "%s\nerror: Step probability must be a "
							"positive constant (line %d).", msg, ln);
				UmbrellaStep *s;
				if (strcmp (tokens[4], "barostat") == 0) {
					Quantity *press = new Quantity (tokens[5], lmp, false, false);
					if (!press->is_valid())
						sprintf (msg, "%s\nerror: \"%s\" is not a valid "
								"pressure (line %d).", msg, tokens[5], ln);
					s = new BarostatStep (lmp, *d, temp, tokens[2], press);
				} else
					s = new UmbrellaStep (lmp, *d, tokens[2]);
				steps_map[tokens[2]] = s;

				if (strcmp (msg, ""))
					return 1;
			} else if (strcmp (tokens[1], "parameter") == 0) {
				if (n != 5)
					sprintf (msg, "Incorrect number of tokens to define "
							"parameter (expected five, line %d).\n", ln);
				Quantity *param = new Quantity (tokens[2], lmp, false,false);
				if (!param->is_valid())
					sprintf (msg, "%s\"%s\" is not a valid parameter name "
							"(line %d).\n", msg, tokens[2], ln);
				Quantity *target = new Quantity (tokens[3],lmp,false,false);
				if (!target->is_valid())
					sprintf (msg, "%s\"%s\" is not a valid target quantity "
							"(line %d).\n", msg, tokens[3], ln);
				// repulsive potentials should work, but let's be safe
				Quantity *spring = new Quantity (tokens[4],lmp, true, false);
				if (!spring->is_valid())
					sprintf (msg, "%s\"%s\" is not a valid spring quantity "
							"(line %d).\n", msg, tokens[4], ln);

				if (strcmp(msg, ""))
					return 1;

				UmbrellaParameter *p = new UmbrellaParameter \
									   (param, target, spring, temp, lmp);
				params.push_back (p);

			} else if (strcmp (tokens[1], "bias_every") == 0) {
				bias_every = (int) std::strtol(tokens[2], &e, 0);
				if (*e != 0) {
					sprintf (msg, "Not an integer! (line %d)\n", ln);
					return 1;
				}
			} else if (strcmp (tokens[1], "take_step") == 0) {
				if (steps_map.find(tokens[2]) == steps_map.end()) {
					sprintf (msg, "Parse error: take_step before %s defined "
							"(line %d).\n", tokens[2], ln);
					return 1;
				}
				current_block = steps_map[tokens[2]]->get_take_step_block();
			} else if (strcmp (tokens[1], "if_accept") == 0) {
				if (steps_map.find(tokens[2]) == steps_map.end()) {
					sprintf (msg, "Parse error: if_accept before %s defined "
							"(line %d).\n", tokens[2], ln);
					return 1;
				}
				current_block = steps_map[tokens[2]]->get_if_accept_block();
			} else if (strcmp (tokens[1], "if_reject") == 0) {
				if (steps_map.find(tokens[2]) == steps_map.end()) {
					sprintf (msg, "Parse error: if_reject before %s defined "
							"(line %d).\n", tokens[2], ln);
					return 1;
				}
				current_block = steps_map[tokens[2]]->get_if_reject_block();
			} else if (strcmp (tokens[1], "step_init") == 0) {
				if (steps_map.find(tokens[2]) == steps_map.end()) {
					sprintf (msg, "Parse error: step_init before %s defined "
							"(line %d).\n", tokens[2], ln);
					return 1;
				}
				current_block = steps_map[tokens[2]]->get_step_init_block();
			} else if (strcmp (tokens[1], "do_every") == 0) {
				int p = (int) std::strtol(tokens[2], &e, 0);
				if (*e != 0) {
					fprintf (stderr, "Not an integer! (line %d)\n", ln);
					return 1;
				}
				PeriodicTask *pt = new PeriodicTask (lmp, p);
				current_block = pt->get_task_block();
				tasks.push_back(pt);
				//delete pt; TODO kill this somewhere
				// Special tasks which get intercepted before LAMMPS
			} else if (strcmp (tokens[1], "histogram") == 0) {
				UmbrellaParameter *p = NULL;
				for (std::vector<UmbrellaParameter *>::iterator it = \
						params.begin(); it != params.end(); ++it)
					if (strcmp (tokens[2], (*it)->get_name()) == 0)
						p = *it;
				if (p == NULL) {
					sprintf (msg, "Unable to locate parameter named \'%s\' "
							"for histogram! (line %d)", tokens[2], ln);
					return 1;
				}
				float min = std::strtod(tokens[3], &e);
				float max = std::strtod(tokens[4], &e);
				int num = (int) std::strtoul(tokens[5], &e, 0);
				int period = (int) std::strtoul(tokens[6], &e, 0);
				if (*e != 0) {
					sprintf (msg, "Error parsing histogram options! "
							"(line %d)", ln);
					return 1;
				}
				Histogram *h = new Histogram (num, min, max, period, p);
				h->set_filename(tokens[7]); //TODO make sure this exists
				histograms.push_back(h);
				//delete h; TODO kill this somewhere


			} else if (strcmp (tokens[1], "get_positions") == 0) {
				strcpy (line, "GET_ATOMS");
			} else if (strcmp (tokens[1], "put_positions") == 0) {
				strcpy (line, "PUT_ATOMS");
			} else if (strcmp (tokens[1], "get_types") == 0) {
				strcpy (line, "GET_TYPES");
			} else if (strcmp (tokens[1], "put_types") == 0) {
				strcpy (line, "PUT_TYPES");
			} else if (strcmp (tokens[1], "force_accept") == 0) {
				strcpy (line, "FORCE_ACCEPT");
				//} else if (strcmp (tokens[1], "do_step") == 0) {
				//sprintf (line, "DO_STEP %s", tokens[2]);  // TODO check third null
		} else {
			sprintf (msg, "Directive not recognized: %s\n", line);
			return 1;
		}
		}
		if (line[0] == '#' || line[0] == '\0') continue; //perhaps pass to LAMMPS so they show up on logs
		current_block->push_back (line);
	}
	delete [] file_data;


	nsteps = steps_map.size();
	steps = new UmbrellaStep *[nsteps];
	float sum = 0;
	int i = 0;
	for (std::map<std::string, UmbrellaStep*>::iterator \
			it = steps_map.begin(); it != steps_map.end(); ++it) {
		it->second->rand_min = sum;
		sum += it->second->probability;
		it->second->rand_max = sum;
		steps[i] = it->second;
		++i;

	}
	if (nsteps && sum != 1) {
		sprintf (msg, "Sum of step probabilities is not one.");
		return 1;
	}

	nparams = params.size();
	param_ptrs = new UmbrellaParameter *[nparams];
	for (int i = 0; i < nparams; ++i)
		param_ptrs[i] = params[i];
	return 0;
}

char *Parser::error_message() {
	return msg;
}

std::string Parser::process_brackets(char *line) {
	//TODO pass line number for error messages?
	char msg[500];
	char file_line[MAX_LINE_LENGTH];
	char file_data[MAX_LINE_LENGTH*Global::get_instance()->get_num_windows()];
	int n = strlen (line);
	char *left = 0;
	char *right = 0;
	int i, j;
	char window_str[100];  //TODO I suppose it could be overrun
	sprintf (window_str, "%d", Global::get_instance()->get_window_index());
	int window_len = strlen (window_str);
	//
	// Substitute @ for rank
	// TODO maybe do this after brackets
	for (i = 0; i < n; ++i)
		if (line[i] == '@') {
			if (window_len > 1)
				for (j = n; j > i; --j)
					line[j + window_len - 1] = line[j];
			// TODO error if line gets overrun
			for (j = 0; j < window_len; ++j)
				line[i + j] = window_str[j];
		}

	// Search for <<
	for (i = 0; i < n - 1; ++i)
		if (line[i] == '<' && line[i + 1] == '<') {
			left = &line[i + 2];
			break;
		}
	// Search for >>
	if (left) {
		for (; i < n - 1; ++i)
			if (line[i] == '>' && line[i + 1] == '>') {
				right = &line[i];
				break;
			}
		if (!right) Global::get_instance()->abort ((char*) \
				"No closing brackets detected!");
		*(left - 2) = *right = '\0';
	} else
		return std::string(line);
	if (right == left) Global::get_instance()->abort ((char*) \
			"Empty brackets encountered in script.");

	if (Global::get_instance()->get_global_rank() == 0) {
		FILE *fp = fopen (left, "r");
		if (fp == NULL) {
			sprintf (msg, "Error opening bracketed file: %s", left);
			Global::get_instance()->abort (msg);
		}
		int i = 0;
		while (!feof (fp)) {
			fgets (file_line, MAX_LINE_LENGTH, fp);
			if (feof (fp)) break;
			if (file_line[0] == '#' || file_line[0] == '\n') continue; //TODO whitespace
			for (unsigned j = 0; j < strlen (file_line); ++j) // Snip newline
				if (file_line [j] == '\n' || file_line [j] == '\r')
					file_line[j] = '\0';
			strcpy (file_data + i * MAX_LINE_LENGTH, file_line);
			++i;
			if (i > Global::get_instance()->get_num_windows()) {
				sprintf (msg, "Number of lines in bracketed file \"%s\""
						" is greater than the number of defined windows (%d). "
						" The first %d lines will be distributed to windows.",\
						left, Global::get_instance()->get_num_windows(), \
						Global::get_instance()->get_num_windows());
				Global::get_instance()->warn(msg);
				break;
			}
		}
		fclose(fp);
		if (i < Global::get_instance()->get_num_windows()) {
			sprintf (msg, "Number of lines in bracketed file \"%s\""
					" is less than the number of defined windows (%d)", \
					left, Global::get_instance()->get_num_windows());
			Global::get_instance()->abort (msg);
		}
	}
	if (Global::get_instance()->get_local_rank() == 0)
		MPI_Scatter (file_data, MAX_LINE_LENGTH, MPI_CHAR, \
				file_line, MAX_LINE_LENGTH, MPI_CHAR, 0, \
				Global::get_instance()->roots_comm);
	MPI_Bcast (file_line, MAX_LINE_LENGTH, MPI_CHAR, 0, \
			Global::get_instance()->local_comm);
	return std::string(line) + std::string(file_line) + std::string(right + 2);
}

void Parser::execute_init() {
	UmbrellaStep::execute_block(lmp, init_block);
}

void Parser::print() {
	fprintf (stdout, "\nInitialization:\n");
	for (unsigned i = 0; i < init_block.size(); ++i)
		fprintf (stdout, "\t%s\n", init_block[i].c_str());
	fprintf (stdout, "Defined steps:\n");
	std::vector<std::string> *v;
	for (std::map<std::string, UmbrellaStep*>::iterator \
			it = steps_map.begin(); it != steps_map.end(); ++it) {
		fprintf ( stdout, "\tStep type: %s\n", it->first.c_str());
		v = it->second->get_take_step_block();
		fprintf (stdout, "\t\tTake step:\n");
		for (unsigned j = 0; j < v->size(); ++j)
			fprintf (stdout, "\t\t\t%s\n", v->at(j).c_str());

		v = it->second->get_if_reject_block();
		fprintf (stdout, "\t\tIf reject:\n");
		for (unsigned j = 0; j < v->size(); ++j)
			fprintf (stdout, "\t\t\t%s\n", v->at(j).c_str());

		v = it->second->get_if_accept_block();
		fprintf (stdout, "\t\tIf accept:\n");
		for (unsigned j = 0; j < v->size(); ++j)
			fprintf (stdout, "\t\t\t%s\n", v->at(j).c_str());

		v = it->second->get_step_init_block();
		fprintf (stdout, "\t\tStep init:\n");
		for (unsigned j = 0; j < v->size(); ++j)
			fprintf (stdout, "\t\t\t%s\n", v->at(j).c_str());
	}
	fprintf (stdout, "Defined parameters:\n");
	for (unsigned j = 0; j < params.size(); ++j)
		fprintf (stdout, "\t%s\n", params[j]->get_name());
}
