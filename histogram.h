#ifndef HISTOGRAMMER_H
#define HISTOGRAMMER_H

#include <stdio.h>

class Histogram {
	public:
		Histogram (int nbins, double min, double max);
		~Histogram ();
		void reset ();
		int update (double val);
		void write (FILE *f);
	private:
		int nbins;
		double min;
		double max;
		unsigned *hist;
};

#endif
