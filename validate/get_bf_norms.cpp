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
//typedef Basis<T, Primal, Interval, Dijkema> IntervalBasis;
typedef Basis<T, Primal, Periodic, CDF> IntervalBasis;

int main(int argc, char *argv[]){

	if(argc != 6){
		cerr << "Usage: " << argv[0] << " d d_ j0 j bc" << endl;
		exit(1);
	}

	int d = atoi(argv[1]);
	int d_ = atoi(argv[2]);
	int j0 = atoi(argv[3]);
	int j = atoi(argv[4]);
	int bc = atoi(argv[5]);

	IntervalBasis basis(d,d_,j0);
	if(bc){
	//	basis.enforceBoundaryCondition<DirichletBC>();

	}

	Integral<Gauss, IntervalBasis, IntervalBasis> integral(basis, basis);
	integral.quadrature.setOrder(20);

	cout << " LEVEL j = " << j  << "               L2-Norm^2 - H1-Norm^2" << endl << endl;

	cout << "------- BSplines -----------" << endl;
	cout << "Left: " << basis.mra.rangeIL(j) << endl;
	for(int k = basis.mra.rangeIL(j).firstIndex(); k <= basis.mra.rangeIL(j).lastIndex(); ++k){
		cout << "k = " << k << ": " <<  integral(j, k, XBSpline, 0, j, k, XBSpline, 0)
			 << " " << integral(j, k, XBSpline, 1, j, k, XBSpline, 1) << endl;
	}

	cout << "Inner: " << basis.mra.rangeII(j) << endl;
	for(int k = basis.mra.rangeII(j).firstIndex(); k <= basis.mra.rangeII(j).lastIndex(); ++k){
		cout << "k = " << k << ": " <<  integral(j, k, XBSpline, 0, j, k, XBSpline, 0)
			 << " " << integral(j, k, XBSpline, 1, j, k, XBSpline, 1) << endl;
	}

	cout << "Right: " << basis.mra.rangeIR(j) << endl;
	for(int k = basis.mra.rangeIR(j).firstIndex(); k <= basis.mra.rangeIR(j).lastIndex(); ++k){
		cout << "k = " << k << ": " <<  integral(j, k, XBSpline, 0, j, k, XBSpline, 0)
			<< " " << integral(j, k, XBSpline, 1, j, k, XBSpline, 1) << endl;
	}

	cout << "------- Wavelets -----------" << endl;
	cout << "Left: " << basis.rangeJL(j) << endl;
	for(int k = basis.rangeJL(j).firstIndex(); k <= basis.rangeJL(j).lastIndex(); ++k){
		cout << "k = " << k << ": " <<  integral(j, k, XWavelet, 0, j, k, XWavelet, 0)
			 << " " << integral(j, k, XWavelet, 1, j, k, XWavelet, 1) << endl;
	}

	cout << "Inner: " << basis.rangeJI(j) << endl;
	for(int k = basis.rangeJI(j).firstIndex(); k <= basis.rangeJI(j).lastIndex(); ++k){
		cout << "k = " << k << ": " <<  integral(j, k, XWavelet, 0, j, k, XWavelet, 0)
			 << " " << integral(j, k, XWavelet, 1, j, k, XWavelet, 1) << endl;
	}

	cout << "Right: " << basis.rangeJR(j) << endl;
	for(int k = basis.rangeJR(j).firstIndex(); k <= basis.rangeJR(j).lastIndex(); ++k){
		cout << "k = " << k << ": " <<  integral(j, k, XWavelet, 0, j, k, XWavelet, 0)
			<< " " << integral(j, k, XWavelet, 1, j, k, XWavelet, 1) << endl;
	}

	return 0;
}


