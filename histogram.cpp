#include <stdio.h>
#include <math.h>

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
	sum = 0;
	sum_squares = 0;
	num_samples = 0;
	for (int i = 0; i < nbins; ++i)
		hist[i] = 0;
}

int Histogram::update (double val) {
	sum += val;
	sum_squares += val*val;
	++num_samples;
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

double Histogram::get_mean() {
	return sum / num_samples;
}

double Histogram::get_standard_deviation() {
	double mean = get_mean();
	return sqrt(sum_squares / num_samples - mean*mean);
}

void Histogram::write (FILE *f) {
	// Write index|count|center
	for (int i = 0; i < nbins; ++i)
		if (hist[i])
			fprintf (f, "%d\t%d\t%f\n", i, hist[i], (min + (i + 0.5) * bin_width) / 3.4); //XXX dirty
			//fprintf (f, "%d\t%d\t%f\n", i, hist[i], min + (i + 0.5) * bin_width);
	fclose (f);
}
