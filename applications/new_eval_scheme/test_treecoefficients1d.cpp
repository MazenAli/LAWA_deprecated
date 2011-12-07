#include <iostream>
#include <lawa/lawa.h>
#include <applications/new_eval_scheme/treecoefficients1d.h>

using namespace std;
using namespace lawa;

typedef double T;

typedef Basis<T,Primal,Interval,Dijkema>                            PrimalBasis;

int main () {

    PrimalBasis basis(2,2,2);

    Coefficients<Lexicographical,T,Index1D> coeff;

    for (int k=basis.mra.rangeI(basis.j0).firstIndex(); k<=basis.mra.rangeI(basis.j0).lastIndex(); ++k) {
        coeff[Index1D(basis.j0,k,XBSpline)] = 1.;
    }
    for (int j=basis.j0; j<=basis.j0+4; ++j) {
        for (int k=basis.rangeJ(j).firstIndex(); k<=basis.rangeJ(j).lastIndex(); ++k) {
            coeff[Index1D(j,k,XWavelet)] = 1.2;
        }
    }

    TreeCoefficients1D<T> treecoefficients;
    treecoefficients = coeff;
    cout << treecoefficients << endl;

    cout << treecoefficients(3) << endl;


    cout << "No segmentation fault so far..." << endl;

    return 0;
}
