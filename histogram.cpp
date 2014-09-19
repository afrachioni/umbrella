#include "histogram.h"

Histogram::Histogram (int nbins, double min, double max) {
	this->nbins = nbins;
	this->min = min;
	this->max = max;
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
	++hist[(int) ((val - min) / max)];
	return 0;
}

void write (FILE *f) {
}
