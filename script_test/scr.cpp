#include <stdio.h>
#include <string.h>

#include <vector>
#include <map>

#include "step_description.h"


int main (int nargs, char **args) {
	FILE *in = fopen ("in.txt", "r");
	if (in == NULL) {
		fprintf (stdout, "Cannot open file!\n");
		return 1;
	}


	int directives = 0;
	int lines = 0;
	char line[500];
	std::vector<char*> *current_block = NULL;
	std::vector<char*> init_block;
	while (1) {
		fgets (line, 500, in);
		if (feof (in)) break;
		if (line[0] != '#' && line[0] != '\n')
			++lines;
		if (strncmp (line, "#AF", 3) == 0)
			++directives;
			

			StepDescription sd ();

		//if (strncmp (line, "#AFINIT", 6) == 0)
			//current_block = &init_block;
		//else if (strncmp (line, "#AFB2", 5) == 0)
			//current_
				

//steptype box_change
//takestep box_change
//ifaccept box_change
//ifreject box_change

	}


	//rewind (in);
	//int i = 0;
	//char lmp_cmd[lines][500];
	//while (1) {
		//fgets (line, 500, in);
		//if (feof (in)) break;
		//if (strncmp (line, "#AF", 3) == 0) {
			//fprintf (stdout, "Special directive: %s", line);
		//}
		//if (line[0] == '#' || line[0] == '\n') continue;
		//strcpy (lmp_cmd[i], line);
		//++i;
	//}

	//for (i = 0; i < lines; ++i)
		//fprintf (stdout, lmp_cmd[i]);


}
