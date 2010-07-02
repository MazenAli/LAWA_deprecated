#include <iostream>
#include <lawa/lawa.h>
#include <extensions/flens/flens.h>

using namespace std;
using namespace lawa;

typedef double T;
typedef SparseGeMatrix<CRS<T,CRS_General> > SparseMatrixT;
typedef flens::DenseVector<flens::Array<T> > DenseVectorT;

int main(){

	DenseVectorT x(4);
	x = 1,1,1,1;
	
	SparseMatrixT A(4,4);
	A(1,1) = 1;
	A(2,2) = 1;
	A(3,4) = 1;
    A(4,1) = 1;
	A.finalize();
	
	DenseVectorT y;
	y = A*x;
	cout << y << endl;

}