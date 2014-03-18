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
typedef Basis<T, Primal, Periodic, CDF> PeriodicBasis;

int main(int argc, char *argv[]){

	if(argc != 5){
		cerr << "Usage: " << argv[0] << " d d_ j0 j" << endl;
		exit(1);
	}

	int d = atoi(argv[1]);
	int d_ = atoi(argv[2]);
	int j0 = atoi(argv[3]);
	int j = atoi(argv[4]);

	PeriodicBasis basis(d,d_,j0);
	
	cout << "INFO:" << endl;
    cout << "  d = " << d << ", d_ = " << d_ << endl;
    cout << "  j0 = " << j0 <<  endl;

	cout << "------- BSplines -----------" << endl;
	cout << j << " " << basis.mra.rangeI(j) << endl;
	cout << "First: " << basis.mra.phi.support(j, basis.mra.rangeI(j).firstIndex()) << endl;
	cout << "Last: " << basis.mra.phi.support(j, basis.mra.rangeI(j).lastIndex()) << endl << endl;

	cout << "Left Indices: " << basis.mra.cardIL(j) << " " << basis.mra.rangeIL(j) << endl;
	cout << "Inner Indices: " << basis.mra.cardII(j) << " " << basis.mra.rangeII(j) << endl;
	cout << "Right Indices: " << basis.mra.cardIR(j) << " " << basis.mra.rangeIR(j) << endl;

    stringstream bsplinefilename;
    bsplinefilename << "bsplines_periodic_d" << d << d_ << "_j0" << j0 << "_" << j << ".txt";  
	ofstream bsplinefile(bsplinefilename.str().c_str());
	for(T x = 0; x <= 1; x += pow2i<T>(-4-j)){
		bsplinefile << x << " ";
		for(int k = basis.mra.rangeI(j).firstIndex(); k <= basis.mra.rangeI(j).lastIndex(); ++k){
			bsplinefile << basis.generator(XBSpline).operator()(x,j,k, 0) << " ";
		}
		bsplinefile << endl;
	}
	bsplinefile.close();

	cout << "------- Wavelets -----------" << endl;
	cout << j << " " << basis.rangeJ(j) << endl;
	cout << "First: " << basis.psi.support(j, basis.rangeJ(j).firstIndex()) << endl;
	cout << "Last: " << basis.psi.support(j, basis.rangeJ(j).lastIndex()) << endl << endl;

	cout << "Left Indices: " << basis.cardJL(j) << " " << basis.rangeJL(j) << endl;
	cout << "Inner Indices: " << basis.cardJI(j) << " " << basis.rangeJI(j) << endl;
	cout << "Right Indices: " << basis.cardJR(j) << " " << basis.rangeJR(j) << endl;


    stringstream waveletfilename;
    waveletfilename << "wavelets_periodic_d" << d << d_ << "_j0" << j0 << "_" << j << ".txt";  
	ofstream waveletfile(waveletfilename.str().c_str());
	for(T x = 0; x <= 1; x += pow2i<T>(-4-j)){
		waveletfile << x << " ";
		for(int k = basis.rangeJ(j).firstIndex(); k <= basis.rangeJ(j).lastIndex(); ++k){
			waveletfile << basis.generator(XWavelet).operator()(x,j,k, 0) << " ";
		}
		waveletfile << endl;
	}




	return 0;
}


