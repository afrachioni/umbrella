#ifndef HISTOGRAM_H
#define HISTOGRAM_H

#include "umbrella_parameter.h"

class Histogram {
	public:
		Histogram (int nbins, double min, double max, unsigned period);
		Histogram (int nbins, double min, double max, unsigned period, \
				UmbrellaParameter *p);
		~Histogram ();
		void reset ();
		int update (double val);
		int update ();
		void write (FILE *f);
		void write (unsigned step_index);
		void write_stats (FILE *f);
		double get_mean(), get_standard_deviation();
		void set_filename (char *filename);
	private:
		void init(int nbins, double min, double max, unsigned period);
		int nbins;
		double min;
		double max;
		UmbrellaParameter *p;
		double bin_width;
		double sum, sum_squares;
		unsigned *hist, num_samples, period;
		char *filename;
};

#endif
