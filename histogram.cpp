#include <stdio.h>

#include "histogram.h"

#include "umbrella_parameter.h"

Histogram::Histogram (int nbins, double min, double max) {
	init (nbins, min, max);
}

Histogram::Histogram (int nbins, double min, double max, UmbrellaParameter *p) {
	init (nbins, min, max);
	this->p = p;
}

void Histogram::init(int nbins, double min, double max) {
	this->nbins = nbins;
	this->min = min;
	this->max = max;
	this->bin_width = (max - min) / nbins;
	hist = new unsigned[nbins];
	reset();
}

Histogram::~Histogram () {
	delete [] hist;
}

void Histogram::reset () {
	for (int i = 0; i < nbins; ++i)
		hist[i] = 0;
}

int Histogram::update (double val) {
	if (val < min)
		return -1;
	else if (val > max)
		return -1;
	++hist[(int) ((val - min) / bin_width)];
	return 0;
}

int Histogram::update() {
	// TODO error if p == NULL
	return update (p->get_last_accepted_value());
}

void Histogram::write (FILE *f) {
	// Write index|count|center
	for (int i = 0; i < nbins; ++i)
		if (hist[i])
			fprintf (f, "%d\t%d\t%f\n", i, hist[i], min + (i - 0.5) * bin_width);
	fclose (f);
}
