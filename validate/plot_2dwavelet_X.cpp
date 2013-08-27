/*
 * plot_2dwavelet_X.cpp
 *
 *  Created on: 14.08.2013
 *      Author: ksteih
 */

#include <iostream>
#include <cstdlib>
#include <lawa/lawa.h>

using namespace std;
using namespace lawa;

typedef double T;
typedef Basis<T, Primal, Periodic, CDF> 	PeriodicBasis;
typedef Basis<T,Orthogonal,Interval,Multi> IntervalBasis;

typedef TensorBasis2D<Adaptive,PeriodicBasis,IntervalBasis>  Basis2D;




int main(int argc, char *argv[]){

	if(argc != 11){
		cerr << "Usage: " << argv[0] << " d d_ j0 outfile Xtype_t(0/1) j_t k_t Xtype_x(0/1) j_x k_x " << endl;
		exit(1);
	}

	int d = atoi(argv[1]);
	int d_ = atoi(argv[2]);
	int j0 = atoi(argv[3]);
	XType xtype_t = atoi(argv[5])==0 ? XBSpline : XWavelet;
	int j_t = atoi(argv[6]);
	int k_t = atoi(argv[7]);
	XType xtype_x = atoi(argv[8])==0 ? XBSpline : XWavelet;
	int j_x = atoi(argv[9]);
	int k_x = atoi(argv[10]);

	PeriodicBasis basis_time(d,d_,j0);
	IntervalBasis basis_space(d,0);

	Basis2D basis2d(basis_time, basis_space);

    std::stringstream outfilename;
    outfilename << argv[4] << "_" << j_t << "_" << k_t << "_" << j_x << "_" << k_x << ".txt";

    std::ofstream outfile(outfilename.str().c_str());
    outfile.precision(16);

    T h = 0.01;

    for (T x=0; x<=1; x+=h) {
        for (T y=0; y<=1; y+=h) {
        	outfile << x << " " << y  << " " << basis2d.first.generator(xtype_t)(x,j_t,k_t,0) * basis2d.second.generator(xtype_x)(y,j_x,k_x,0)  << std::endl;
        }
        outfile << std::endl;
    }
    outfile.close();

	return 0;
}







