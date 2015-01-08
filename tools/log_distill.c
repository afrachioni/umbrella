#include<stdio.h>

// Usage: log_distill <in_file> <out_file> <stride>

int main (int nargs, char **args) {
	FILE *in = fopen (args[1], "r");
	FILE *out = fopen (args[2], "w");
	int count = 0;
	int stride = atoi(args[3]);
	char line[200];
	fgets (line, 200, in);

	while (!feof (in)) {
		if (count % stride == 0 || line[0] == '#')
			fprintf (out, line);
		fgets (line, 200, in);
		++count;
	}
	fclose(in); fclose(out);
	return 0;
}
