#include <stdio.h>
#include <string.h>

#include <vector>
#include <map>
#include <string>

#include "step_description.h"


int main (int nargs, char **args) {
	FILE *in = fopen ("in.txt", "r");
	if (in == NULL) {
		fprintf (stderr, "Cannot open file!\n");
		return 1;
	}


//steptype box_change
//takestep box_change
//ifaccept box_change
//ifreject box_change

	std::map<std::string, StepDescription> steps_map;
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
	while (1) {
		++lnum;
		fgets (line, 500, in);
		if (feof (in)) break;
		//fprintf (stdout, line);
		n = sscanf (line, "%s %s %s %s", first_token, second_token, third_token, fourth_token);
		if (strcmp (first_token, "#AF") == 0 && n > 0) {
			fprintf (stdout, "Special directive: %s", line);
			if (n == 1) {
				fprintf (stderr, "Parse error: empty directive at line %d.\n", lnum);
				break;
			} else if (strcmp (second_token, "steptype") == 0) {
				fprintf (stdout, "Type: %s\n", third_token);
				steps_map[third_token];

			} else if (strcmp (second_token, "takestep") == 0) {

				fprintf (stdout, "Takestep: %s\n", third_token);
				if (steps_map.find(third_token) == steps_map.end()) {
					fprintf (stderr, "Parse error: takstep before %s defined\n", third_token);
					break;
				}
				fprintf (stdout, "Setting current block to take step for %s\n", third_token);
				current_block = steps_map[third_token].get_take_step_block();
			} else {
				fprintf (stdout, "Directive not recognized: %s\n", line);
			}
		}
		if (line[0] == '#' || line[0] == '\n') continue; //perhaps pass these to LAMMPS so they show up on the logs

		if (current_block != NULL) {
			current_block->push_back (line);
			fprintf (stdout, "pushed back: %s", line);
		}
	}

	for (std::map<std::string, StepDescription>::iterator it = steps_map.begin(); it != steps_map.end(); ++it) {
		fprintf ( stdout, "Map key: %s\n", it->first.c_str());
		std::vector<std::string> *v = it->second.get_take_step_block();
		fprintf (stdout, "num take_step lines: %d\n", v->size());
		for (int j = 0; j < v->size(); ++j)
			fprintf (stdout, "\tTake step values: %s", (*v)[j].c_str());//v->at(j));
	}
}
