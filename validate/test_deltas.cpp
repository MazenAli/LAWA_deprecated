#include <iostream>
#include <fstream>
#include <lawa/lawa.h>

using namespace std;
using namespace lawa;

typedef double T;
typedef flens::GeMatrix<flens::FullStorage<T, cxxblas::ColMajor> >  FullColMatrixT;

typedef Basis<T, Primal, Interval, Dijkema>                         PrimalBasis;

int main()
{
    /// wavelet basis parameters:
    int d = 2;          // (d,d_)-wavelets
    int d_ = 2;
    int j0 = 2;         // minimal level
    int J = 5;          // maximal level

    PrimalBasis basis(d, d_, j0);
    basis.enforceBoundaryCondition<DirichletBC>();
    int j = 3, k = 4;

    FullColMatrixT deltas = computeDeltas<T,PrimalBasis>(basis, j, k, XBSpline);
    cout << deltas << endl;

    return 0;

}
