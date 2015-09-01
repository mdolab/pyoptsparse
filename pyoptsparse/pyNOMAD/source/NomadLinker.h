#ifndef NOMADLINKER_H
#define NOMADLINKER_H

#include <vector>

class Callback;

class NomadLinker
{
    public:
        NomadLinker();
        ~NomadLinker();
        void setCallback( Callback &callback);
		std::vector<double> call( int n, int m, double* x, int xdim, double* lb, int lbdim,
				double* ub, int ubdim, double min_poll_size, double min_mesh_size, int max_bbe, int display_degree, int print_file );
        void setSolution( std::vector<double> transfer);
		std::vector<double> getSolution();
		Callback* passcallback() { return callback_; }

    private:
		std::vector<double> solution;
        Callback* callback_;
};

#endif
