#include <iostream>
#include <lawa/lawa.h>

using namespace std;
using namespace lawa;

typedef double T;

typedef Basis<T,Primal,Interval,SparseMulti> Basis1D;

typedef Integral<Gauss,Basis1D,Basis1D>     Integral1D;

int main (int argc, char *argv[]) {

    cout.precision(16);
    int d=4;
    int j0=2;

    Basis1D basis(d,j0);
    basis.enforceBoundaryCondition<DirichletBC>();

    int j=2;

    Integral1D integral(basis,basis);



    cout << "Basis: range = " << basis.rangeJ(j) << endl;
    for(int k=basis.rangeJ(j).firstIndex(); k<=basis.rangeJ(j).lastIndex(); ++k) {
        ofstream file("sparsemulti_wavelet.dat");
        cout << basis.psi.support(j,k) << endl;
        cout << basis.psi.singularSupport(j,k) << endl;
        cout << "Integral-test: " << integral(j,k,XWavelet,0,j,k,XWavelet,0) << endl;
        for (T x=-2.; x<=2.; x+=pow2i<T>(-9)) {
            T h = 1e-4;
            T approx_d_psi = (basis.psi(x+h,j,k,0)-basis.psi(x-h,j,k,0))/(2.*h);
            file << x << " " << basis.psi(x,j,k,0) << " " << basis.psi(x,j,k,1) << " " << approx_d_psi << endl;
        }
        file.close();
        getchar();
    }

    cout << "MRA: range = " << basis.mra.rangeI(j) << endl;
    for (int k=basis.mra.rangeI(j).firstIndex(); k<=basis.mra.rangeI(j).lastIndex(); ++k) {
        ofstream file("sparsemulti_scaling.dat");

        cout << basis.mra.phi.support(j,k) << endl;
        cout << basis.mra.phi.singularSupport(j,k) << endl;
        cout << "Integral-test: " << integral(j,k,XBSpline,0,j,k,XBSpline,0) << endl;
        for (T x=-2.; x<=2.; x+=pow2i<T>(-8)) {
            T h = 1e-4;
            T approx_d_phi = (basis.mra.phi(x+h,j,k,0)-basis.mra.phi(x-h,j,k,0))/(2.*h);
            file << x << " " << basis.mra.phi(x,j,k,0) << " " << basis.mra.phi(x,j,k,1) << " " << approx_d_phi << endl;
        }
        file.close();
        getchar();
    }

    return 0;
}
