/*
 * test_index.cc
 *
 *  Created on: 09.08.2010
 *      Author: sebastian
 */

#include <iostream>
#include <adaptive/index.h>
#include <lawa/lawa.h>

typedef double T;

using namespace lawa;
using namespace std;

int main()
{
	int d=2,d_=2;

	Index index1(d,d_);
	Index index2(d,d_,3,4,XWavelet);
	cout << "Index1 = " << index1 << endl;
	index1 = index2;
	cout << "Index1 after assignement: " << index1 << endl;
	Basis<T,Primal,Interval,Dijkema> basis(d,d_,4);
	IntervalIndex<Basis<T,Primal,Interval,Dijkema> > index3(basis,d,d_);
	cout << index3 << endl;

	return 0;


}
