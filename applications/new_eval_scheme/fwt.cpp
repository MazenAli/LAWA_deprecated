#include <iostream>
#include <lawa/lawa.h>

using namespace std;
using namespace lawa;

typedef double T;

typedef Basis<T,Primal,Interval,Dijkema>                            PrimalBasis;
typedef Basis<T,Dual,Interval,Dijkema>                              DualBasis;

typedef flens::DenseVector<flens::Array<T> >                        DenseVectorT;

int main (int argc, char *argv[]) {

    if (argc!=5) {
        cout << "Usage: " << argv[0] << " d d_ j0 J" << endl;
        return 0;
    }
    int d  = atoi(argv[1]);
    int d_ = atoi(argv[2]);
    int j0 = atoi(argv[3]);
    int J  = atoi(argv[4]);

    PrimalBasis basis(d,d_,j0);
    DualBasis   dual_basis(d,d_,j0);

    DenseVectorT c(basis.mra.rangeI(J));
    for (int i=basis.mra.rangeI(j0).firstIndex(); i<=basis.mra.rangeI(j0).lastIndex(); ++i) {
        c(i) = T(rand()) / RAND_MAX;
    }
    for (int j=j0; j<=J-1; ++j) {
        for (int i=basis.rangeJ(j).firstIndex(); i<=basis.rangeJ(j).lastIndex(); ++i) {
            c(basis.mra.cardI(j)+i) = T(rand()) / RAND_MAX;;
        }
    }

    DenseVectorT c_single(basis.mra.rangeI(J));
    ifwt(c, basis, J-1, c_single);

    DenseVectorT c_multi(basis.mra.rangeI(J));
    fwt(c_single, dual_basis, J-1, c_multi);


    DenseVectorT c_diff(basis.mra.rangeI(J));
    c_diff = c - c_multi;
    long double ell2norm = 0.L;
    for (int i=c_diff.firstIndex(); i<=c_diff.lastIndex(); ++i) {
        long double tmp = (long double)c_diff(i);
        ell2norm += tmp*tmp;
    }
    ell2norm = sqrt(ell2norm);
    cout << "|| c_diff || = " << ell2norm << endl;


/*
    for (int l=basis.j0; l<=J; ++l) {
        PrimalBasis basis_shifted(d,d_,l);

        stringstream filename;
        filename << "u_shift_" << l-basis.j0 << ".dat";
        ofstream file(filename.str().c_str());
        for (T x=0.; x<=1.; x+=pow2i<T>(-J-2)) {
            file << x << " " << evaluate(basis_shifted, J, c, x, 0) << endl;
        }
        file.close();

        if (l==J) break;

        DenseVectorTView cview = c(basis.mra.rangeI(l+1));
        DenseVectorT     z = c(basis.mra.rangeI(l+1));
        reconstruct(z, basis, l, cview);
    }
*/


    return 0;
}
