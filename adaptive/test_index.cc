/*
 * test_index.cc
 *
 *  Created on: 09.08.2010
 *      Author: sebastian
 */

#include <iostream>
#include <adaptive/index.h>
#include <adaptive/indexset.h>
#include <lawa/lawa.h>


typedef double T;

using namespace lawa;
using namespace std;


int main()
{

	Index1d index1;
	Index1d index2(3,4,XWavelet);

	lt<Lexicographical,Index1d> compare;
	if (compare(index1,index2))	cout << index1 << endl;
	else						cout << index2 << endl;

	IndexSet<Index1d> indexset1, indexset2;
	indexset1.insert(index1); indexset1.insert(Index1d()); indexset1.insert(index2);
	cout << indexset1 << endl;
	Basis<T,Primal,Interval,Dijkema> basis_interval(2,2,2);
	Basis<T,Primal,Periodic,CDF> basis_periodic(2,2,2);
	indexset2 = C(indexset1,0.5,basis_interval);
	indexset2 = C(indexset1,0.5,basis_periodic);

	return 0;

}
