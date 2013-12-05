#include <stdio.h>
#include <string.h>

#include <vector>
#include <map>
#include <string>

#include "umbrella_step.h"
#include "parser.h"


int main (int nargs, char **args) {
	FILE *in = fopen ("in.txt", "r");
	if (in == NULL) {
		fprintf (stderr, "Cannot open file!\n");
		return 1;
	}

	Parser *p = new Parser(in);
	p->parse();
	p->print();
	delete p;
}

Parser::Parser(FILE *in) {
	this->in = in;
};

Parser::~Parser() {};

void Parser::parse() {
	rewind (in);
	int i = 0;
	int n;
	int lnum = 0;
	char line[500];
	char first_token[20];
	char second_token[20];
	char third_token[20];
	char fourth_token[20];

	std::vector<std::string>* current_block = NULL;
	while (1) { // no exceptions yet for too few tokens or things not specified
		++lnum;
		fgets (line, 500, in);
		if (feof (in)) break;
		n = sscanf (line, "%s %s %s %s", first_token, second_token, third_token, fourth_token);
		if (strcmp (first_token, "#AF") == 0 && n > 0) {
			fprintf (stdout, "Special directive: %s", line);
			if (n == 1) {
				fprintf (stderr, "Parse error: empty directive at line %d.\n", lnum);
				break;
			} else if (strcmp (second_token, "steptype") == 0) {
				fprintf (stdout, "Type: %s\n", third_token);
				steps_map[third_token];

			} else if (strcmp (second_token, "spring_var") == 0) {
				if (steps_map.find (third_token) == steps_map.end()) {
					fprintf (stderr, "Parse error: spring_var before %s defined\n", third_token);
					break;
				}
				steps_map[third_token].set_spring_var(fourth_token);
			
			} else if (strcmp (second_token, "target_var") == 0) {
				if (steps_map.find (third_token) == steps_map.end()) {
					fprintf (stderr, "Parse error: target_var before %s defined\n", third_token);
					break;
				}
				steps_map[third_token].set_target_var(fourth_token);

			} else if (strcmp (second_token, "takestep") == 0) {
				fprintf (stdout, "Takestep: %s\n", third_token);
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
				fprintf (stdout, "Directive not recognized: %s\n", line);
			}
		}
		if (line[0] == '#' || line[0] == '\n') continue; //perhaps pass to LAMMPS so they show up on the logs
		if (current_block != NULL) {
			current_block->push_back (line);
			fprintf (stdout, "pushed back: %s", line);
		}
	}


}
void Parser::print() {
	fprintf (stdout, "\n\n\nDefined steps:\n");

	std::vector<std::string> *v;
	for (std::map<std::string, UmbrellaStep>::iterator it = steps_map.begin(); it != steps_map.end(); ++it) {
		fprintf ( stdout, "Step type: %s\n", it->first.c_str());
		fprintf ( stdout, "Spring variable name: %s\n", it->second.get_spring_var());
		fprintf ( stdout, "Target variable name: %s\n", it->second.get_target_var());
		v = it->second.get_take_step_block();
		fprintf (stdout, "\tTake step:\n");
		for (int j = 0; j < v->size(); ++j)
			fprintf (stdout, "\t\t%s", (*v)[j].c_str());//v->at(j));

		v = it->second.get_if_reject_block();
		fprintf (stdout, "\tIf reject:\n");
		for (int j = 0; j < v->size(); ++j)
			fprintf (stdout, "\t\t%s", v->at(j).c_str());

		v = it->second.get_if_accept_block();
		fprintf (stdout, "\tIf accept:\n");
		for (int j = 0; j < v->size(); ++j)
			fprintf (stdout, "\t\t%s", v->at(j).c_str());

		v = it->second.get_step_init_block();
		fprintf (stdout, "\tStep init:\n");
		for (int j = 0; j < v->size(); ++j)
			fprintf (stdout, "\t\t%s", v->at(j).c_str());
	}
}
