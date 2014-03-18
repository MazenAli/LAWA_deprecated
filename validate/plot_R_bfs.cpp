/*
 * plot_periodic_bfs.cpp
 *
 *  Created on: 26.10.2012
 *      Author: ksteih
 */

#include <iostream>
#include <cstdlib>
#include <lawa/lawa.h>

using namespace std;
using namespace lawa;

typedef double T;
typedef Basis<T, Primal, R, CDF> ReallineBasis;
typedef Basis<T, Dual, R, CDF>   ReallineDualBasis;

int main(int argc, char *argv[]){

	if(argc != 7){
		cerr << "Usage: " << argv[0] << " d d_ j0 j kfirst klast" << endl;
		exit(1);
	}

	int d = atoi(argv[1]);
	int d_ = atoi(argv[2]);
	int j0 = atoi(argv[3]);
	int j = atoi(argv[4]);
	int kfirst = atoi(argv[5]);
	int klast = atoi(argv[6]);

	ReallineBasis basis(d,d_,j0);
    ReallineDualBasis dualbasis(d,d_,j0);
	
	cout << "INFO:" << endl;
    cout << "  d = " << d << ", d_ = " << d_ << endl;
    cout << "  j0 = " << j0 <<  endl;
    
	cout << "------- Generators -----------" << endl;
	stringstream generatorfilename;
    generatorfilename << "generators_R_d" << d << d_ << "_j0" << j0 << "_k" << 0 << ".txt";
	ofstream generatorfile(generatorfilename.str().c_str());
	for(T x = dualbasis.generator(XBSpline).support(j0,0).l1; x <= dualbasis.generator(XBSpline).support(j0,0).l2; x += pow2i<T>(-8-j)){
		generatorfile << x << " ";
			generatorfile << basis.generator(XBSpline).operator()(x,j0,0, 0) << " ";
			generatorfile << basis.generator(XWavelet).operator()(x,j0,0, 0) << " ";
			generatorfile << dualbasis.generator(XBSpline).operator()(x,j0,0, 0) << " ";
			generatorfile << dualbasis.generator(XWavelet).operator()(x,j0,0, 0) << " ";
		generatorfile << endl;
	}
	generatorfile.close();
	
	
	cout << "------- BSplines -----------" << endl;
	cout << "j = " << j << endl;
    cout << "k = " << kfirst << ", ..., " << klast << endl;

    stringstream bsplinefilename;
    bsplinefilename << "bsplines_R_d" << d << d_ << "_j0" << j0 << "_" << j << "_k" << kfirst <<"_"<< klast << ".txt";  
	ofstream bsplinefile(bsplinefilename.str().c_str());
	for(T x = basis.generator(XBSpline).support(j,kfirst).l1; x <= basis.generator(XBSpline).support(j,klast).l2; x += pow2i<T>(-4-j)){
		bsplinefile << x << " ";
		for(int k = kfirst; k <= klast; ++k){
			bsplinefile << basis.generator(XBSpline).operator()(x,j,k, 0) << " ";
		}
		bsplinefile << endl;
	}
	bsplinefile.close();

	cout << "------- Wavelets -----------" << endl;
	cout << "j = " << j << endl;
    cout << "k = " << kfirst << ", ..., " << klast << endl;


    stringstream waveletfilename;
    waveletfilename << "wavelets_R_d" << d << d_ << "_j0" << j0 << "_" << j << "_k" << kfirst <<"_" << klast << ".txt";  
	ofstream waveletfile(waveletfilename.str().c_str());
	for(T x = basis.generator(XWavelet).support(j,kfirst).l1; x <= basis.generator(XWavelet).support(j,klast).l2; x += pow2i<T>(-4-j)){
		waveletfile << x << " ";
		for(int k = kfirst; k <= klast; ++k){
			waveletfile << basis.generator(XWavelet).operator()(x,j,k, 0) << " ";
		}
		waveletfile << endl;
	}

	return 0;
}


