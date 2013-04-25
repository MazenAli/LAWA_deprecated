/*
 * test_sparseTransp.cpp
 *
 *  Created on: 16.04.2013
 *      Author: ksteih
 */

#include <iostream>
#include <algorithm>
#include <lawa/lawa.h>

typedef double T;
typedef flens::SparseGeMatrix<flens::CRS<T,flens::CRS_General> >    SparseMatrixT;
typedef flens::GeMatrix<flens::FullStorage<T, cxxblas::ColMajor> >  FullColMatrixT;

using namespace std;
using namespace lawa;

int main (){

	SparseMatrixT A(4,4);
	A(2,3) = 1.;
	A(4,1) = 2.;
	A(3,3) = 3;
	A(1,2) = 4;
	A(1,4) = 5;
	A.finalize();

	cout << "A:" << endl;
	cout << "rows = " << A.engine().rows;
	cout << "cols = " << A.engine().columns;
	cout << "vals = " << A.engine().values;

	SparseMatrixT ApAT(4,4);

	for(int i = 0; i < 2; ++i){
	for(auto it = A.begin(); it != A.end(); ++it){
		// Find transposed entry A_{c,r} to it = A_{r,c}
		int r = (*it).first.first;
		int c = (*it).first.second;

		cout << "Index: (" << r << ", " << c << ")" << endl;

		// Test all elements in row c
		T val = 0;
		for(int k = it._crs.rows(c); k < it._crs.rows(c+1); ++k){
			if(it._crs.columns(k) == r){
				val = it._crs.values(k);
				break;
			}
		}
		T entry = (*it).second + val;
		cout << "Entry = " << entry << endl;
		if(val > 0){
			entry *= 0.5;
		}
		ApAT((*it).first.first, (*it).first.second) += entry;
		ApAT((*it).first.second, (*it).first.first) += entry;
	}
	}

	ApAT.finalize();
	ApAT *= 0.5;

	FullColMatrixT Adense(4,4);
	densify(cxxblas::NoTrans, ApAT, Adense);
	cout << Adense << endl;



	return 0;
}

