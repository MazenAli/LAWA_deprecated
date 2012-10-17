#include <iostream>
#include <lawa/lawa.h>

using namespace std;
using namespace lawa;

typedef double T;

typedef flens::DenseVector<flens::Array<T> >        DenseVectorT;
typedef flens::SparseGeMatrix<CRS<T,CRS_General> >  SparseMatrixT;

/// Basis definitions
typedef Basis<T,Primal,Interval,Dijkema>            PrimalBasis;


T
u(T x) {
    return 20*x*x*(1-x)*(1-x);
}

T
du(T x) {
    return 40*(x*(1-x)*(1-x) - x*x*(1-x));
}

T
ddu(T x) {
    return 40*((1-x)*(1-x) - 4*x*(1-x) + x*x);
}

int main (int argc, char *argv[]) {

    PrimalBasis basis(3,3,3);
    int J = 4;

    int N = basis.mra.cardI(J);

    DenseVectorT sing_pts;
    Function<T>  fct_u(u,sing_pts);
    IntegralF<Gauss,PrimalBasis> integral_u_psi(fct_u,basis);
    integral_u_psi.quadrature.setOrder(8);

    DenseVectorT rhs(basis.mra.rangeI(J)), coeffs(basis.mra.rangeI(J));
    for (int k=basis.mra.rangeI(J).firstIndex(); k<=basis.mra.rangeI(J).lastIndex(); ++k) {
        rhs(k) = integral_u_psi(J,k,XBSpline,0);
    }

    SparseMatrixT A(N,N);
    int offset = basis.mra.rangeI(J).firstIndex() - 1;
    Integral<Gauss,PrimalBasis,PrimalBasis> integral_psi_psi(basis,basis);
    for (int k1=basis.mra.rangeI(J).firstIndex(); k1<=basis.mra.rangeI(J).lastIndex(); ++k1) {
        for (int k2=basis.mra.rangeI(J).firstIndex(); k2<=basis.mra.rangeI(J).lastIndex(); ++k2) {
            T val = integral_psi_psi(J,k1,XBSpline,0,J,k2,XBSpline,0);
            if (fabs(val)>0) A(k1-offset,k2-offset) = val;
        }
    }
    A.finalize();

    int numOfIters = lawa::cg(A,coeffs,rhs, 1e-8, 100);
    cout << "Number of iterations: " << numOfIters << endl;

    ofstream plotfile("visualize_idea_new_residual.dat");

    for (T x=0.; x<=1.; x+=0.0001) {
        T val1 = 0., val2 = 0.;
        for (int k=basis.mra.rangeI(J).firstIndex(); k<=basis.mra.rangeI(J).lastIndex(); ++k) {
            val1 += coeffs(k) * basis.mra.phi(x,J,k,0);
            val2 += coeffs(k) * basis.mra.phi(x,J,k,2);
        }
        plotfile << x << " " << u(x) << " " << val1 << " " << ddu(x) << " " << val2 << endl;
    }

    return 0;
}
