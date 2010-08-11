/*
 * test_index.cc
 *
 *  Created on: 09.08.2010
 *      Author: sebastian
 */

#include <iostream>
#include <adaptive/index.h>
#include <adaptive/indexset.h>
#include <adaptive/coefficients.h>
#include <lawa/lawa.h>


typedef double T;

using namespace lawa;
using namespace std;

struct Timer
{
    void
    start()	{	_start = clock(); }

    void
    stop() {	_stop = clock(); }

    double
    elapsed() {	return double(_stop - _start) / CLOCKS_PER_SEC; }

    clock_t _start;
    clock_t _stop;
};


int main()
{

	cout.precision(16);
	int d=2,d_=2;
	Index1d index1;
	Index1d index2(3,4,XWavelet);
	Index1d index3(4,5,XWavelet);

	lt<Lexicographical,Index1d> compare;
	if (compare(index1,index2))	cout << index1 << endl;
	else						cout << index2 << endl;

	IndexSet<Index1d> indexset1(d,d_), indexset2(d,d_);
	indexset1.insert(index1); indexset1.insert(Index1d()); indexset1.insert(index2);
	cout << indexset1 << endl;
	Basis<T,Primal,Interval,Dijkema> basis_interval(2,2,2);
	Basis<T,Primal,Periodic,CDF> basis_periodic(2,2,2);

	indexset2 = C_realline(indexset1,0.5);
	cout << "C(indexset1,0.5) = " << indexset2 << endl;

	Coefficients<Lexicographical,T,Index1d> u(d,d_), v(d,d_), w;
	u[index1] = 1.;
	u[index2] = 2.;
	v[index2] = 2.;
	v[index3] = 3.;
	cout << "u = " << u << endl;
	cout << "v = " << v << endl;
	cout << "Norm of u = " << u.norm() << endl;
	v = 2.*v;
	cout << "2*v = " << v << endl;
	w = u +v;
	cout << "u+2*v = " << w << endl;

	FillWithZeros(indexset2, w);
	cout << "w filled with zeros = " << w << endl;

	BSpline<T,Primal,R,CDF> phi(d,0);
	Wavelet<T,Primal,R,CDF> psi(d,d_,0);
	Timer time;
	time.start();
	IndexSet<Index1d> Lambda1 = lambdaTilde1d_PDE(Index1d(0,0,XBSpline), phi, psi, 25, 0, 25, false);
	IndexSet<Index1d> Lambda2 = lambdaTilde1d_PDE(Index1d(0,0,XWavelet), phi, psi, 25, 0, 25, false);
	time.stop();
	T lambdaTildeTime = time.elapsed();

	//cout << "Output of lamdaTilde for BSpline with support " << phi.support(0,0)  << " = " << Lambda << endl;
	cout << "Size of lambdaTilde: " << Lambda1.size() << endl;
	cout << "Size of lambdaTilde: " << Lambda2.size() << endl;
	cout << "Time for lambdaTilde: " << lambdaTildeTime << endl;

	Coefficients<AbsoluteValue,T,Index1d> u_abs(d,d_);
	u_abs = u;
	cout << "u_abs = " << u_abs << endl;


	return 0;

}
