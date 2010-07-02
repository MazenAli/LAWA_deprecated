#include <iostream>
#include <lawa/lawa.h>

using namespace std;
using namespace lawa;

typedef double T;
typedef SparseGeMatrix<CRS<T,CRS_General> > SparseMatrixT;
typedef flens::DenseVector<flens::Array<T> > DenseVectorT;

int main(){

	DenseVectorT x(8);
	x = 1,1,1,1,1,1,1,1;
	
	SparseMatrixT A(8,8);
	A(1,1) = 1;
	A(1,5) = 1;
	A(3,4) = 1;
	A.finalize();
	
	DenseVectorT y;
	y = A*x;
	cout << y << endl;

}