#ifndef HISTOGRAM_H
#define HISTOGRAM_H

#include "umbrella_parameter.h"

class Histogram {
	public:
		Histogram (int nbins, double min, double max);
		Histogram (int nbins, double min, double max, UmbrellaParameter *p);
		~Histogram ();
		void reset ();
		int update (double val);
		int update ();
		void write (FILE *f);
	private:
		void init(int nbins, double min, double max);
		int nbins;
		double min;
		double max;
		UmbrellaParameter *p;
		double bin_width;
		unsigned *hist;
};

#endif
