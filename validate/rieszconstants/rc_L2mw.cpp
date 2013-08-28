/*
 * rc_periodic.cpp
 *
 *  Created on: 28.08.2013
 *      Author: ksteih
 */

#include <iostream>
#include <cstdlib>
#include <lawa/lawa.h>

using namespace std;
using namespace lawa;

typedef double T;
typedef Basis<T, Orthogonal,Interval,Multi> L2MWBasis;



typedef AdaptiveIdentityOperator1D<T,Orthogonal,Interval,Multi>   			IdentityOp;

typedef flens::SparseGeMatrix<flens::CRS<T,flens::CRS_General> >    		SparseMatrixT;
typedef flens::SparseSyMatrix<flens::CRS<T,flens::CRS_UpperTriangular> >    SparseSymMatrixT;
typedef flens::DenseVector<flens::Array<T> >                        		DenseVectorT;
typedef flens::GeMatrix<flens::FullStorage<T, cxxblas::ColMajor> >  		FullColMatrixT;


int main(int argc, char *argv[]){

	if(argc != 4){
		cerr << "Usage: " << argv[0] << " d j0 J" << endl;
		exit(1);
	}

	int d = atoi(argv[1]);
	int j0 = atoi(argv[2]);
	int J = atoi(argv[3]);

	cout << "L2 MW Basis with d = " << d << ", j0 = " << j0 << ", J = " << J << endl << endl;

	L2MWBasis basis(d,j0);
	IdentityOp id_op(basis);


	IndexSet<Index1D> Lambda;
	getFullIndexSet(basis, Lambda, J);
	int N = Lambda.size();

    Integral<Gauss,L2MWBasis,L2MWBasis> _integral(basis, basis);


    cout << "--------- L2 -------------" << endl << endl;
	{
    	// Assemble L2-Matrix
    	cout << "Assembling... " << endl;

		FullColMatrixT I_L2(N,N);
		int col_count = 1;
		for (auto& ind_col : Lambda) {
			int row_count = 1;
			T P_col = 1./std::sqrt(_integral(ind_col.j,ind_col.k,ind_col.xtype,0,ind_col.j,ind_col.k,ind_col.xtype,0));
			for (auto& ind_row : Lambda) {
				I_L2(row_count, col_count) = P_col * id_op(ind_row, ind_col) * 1./std::sqrt(_integral(ind_row.j,ind_row.k,ind_row.xtype,0,ind_row.j,ind_row.k,ind_row.xtype,0));
				row_count++;
			}
			col_count++;
		}
		cout << "... done: " << endl;


		cout << "Finding Eigenvalues... " << endl;
		DenseVectorT    wr(N), wi(N);
		FullColMatrixT  vl, vr;
		int info = ev(false, false, I_L2, wr, wi, vl, vr);
		cout << "EV Info: " << info << endl;

		cout << "Largest Eigenvalue: " << wr(1) << endl;
		cout << "Smallest Eigenvalue: " << wr(N) << endl;

		//cout << "WR:" << wr << endl;
		//cout << "WI:" << wi << endl;
	}
	cout << endl;

    cout << "--------- H1 -------------" << endl << endl;
	{
    	// Assemble L2-Matrix
    	cout << "Assembling... " << endl;

		FullColMatrixT I_H1(N,N);
		int col_count = 1;
		for (auto& ind_col : Lambda) {
			int row_count = 1;
			T P_col = 1./std::sqrt(_integral(ind_col.j,ind_col.k,ind_col.xtype,0,ind_col.j,ind_col.k,ind_col.xtype,0)
						         + _integral(ind_col.j,ind_col.k,ind_col.xtype,1,ind_col.j,ind_col.k,ind_col.xtype,1));
			for (auto& ind_row : Lambda) {
				I_H1(row_count, col_count) = P_col * id_op(ind_row, ind_col) *
						1./std::sqrt(_integral(ind_row.j,ind_row.k,ind_row.xtype,0,ind_row.j,ind_row.k,ind_row.xtype,0)
								   + _integral(ind_row.j,ind_row.k,ind_row.xtype,1,ind_row.j,ind_row.k,ind_row.xtype,1));
				row_count++;
			}
			col_count++;
		}
		cout << "... done: " << endl;


		cout << "Finding Eigenvalues... " << endl;
		DenseVectorT    wr(N), wi(N);
		FullColMatrixT  vl, vr;
		int info = ev(false, false, I_H1, wr, wi, vl, vr);
		cout << "EV Info: " << info << endl;

		cout << "Largest Eigenvalue: " << wr(1) << endl;
		cout << "Smallest Eigenvalue: " << wr(N) << endl;

		//cout << "WR:" << wr << endl;
		//cout << "WI:" << wi << endl;
	}


	return 0;
}

