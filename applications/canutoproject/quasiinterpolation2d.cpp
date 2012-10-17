#include <iostream>
#include <fstream>
#include <lawa/lawa.h>
#include <applications/finance/initialconditions/initialconditions.h>

using namespace std;
using namespace lawa;

typedef double T;

typedef flens::DenseVector<flens::Array<T> >                        DenseVectorT;
typedef flens::GeMatrix<flens::FullStorage<T, cxxblas::ColMajor> >  DenseMatrixT;

typedef Basis<T,Primal,Interval,Dijkema>                            PrimalBasis;

typedef Coefficients<Lexicographical,T,Index2D>::const_iterator     const_coeff2d_it;

T f(T x, T y) { return (x+1)*std::exp(x) * exp(y*y); }

int main (int argc, char *argv[]) {

    PrimalBasis basis(3,3,3);
    int J = 4;
    int N = basis.mra.cardI(J);

    T h = 1./(N-1);
    int offset = basis.mra.rangeI(J).firstIndex();

    cout << "N = " << N << ", offset = " << offset << endl;
/*
    Coefficients<Lexicographical,T,Index2D> u;
    for (int k1=basis.mra.rangeI(J).firstIndex(); k1<=basis.mra.rangeI(J).lastIndex(); ++k1) {
        T x1 = (k1-offset)*h;
        Index1D index1(J,k1,XBSpline);
        for (int k2=basis.mra.rangeI(J).firstIndex(); k2<=basis.mra.rangeI(J).lastIndex(); ++k2) {
            T x2 = (k2-offset)*h;
            Index1D index2(J,k2,XBSpline);
            Index2D index(index1,index2);
            u[index] = f(x1,x2) / (basis.mra.phi(x1,J,k1,0) * basis.mra.phi(x2,J,k2,0));
        }
    }

    ofstream plotfile("interpolation.dat");
    for (T x1=0.; x1<=1.; x1+=h) {
        for (T x2=0.; x2<=1.; x2+=h) {
            T val = 0.;
            for (const_coeff2d_it it=u.begin(); it!=u.end(); ++it) {
                int k1 = (*it).first.index1.k, k2 = (*it).first.index2.k;
                val += (*it).second * basis.mra.phi(x1,J,k1,0) * basis.mra.phi(x2,J,k2,0);
            }
            plotfile << x1 << " " << x2 << " " << f(x1,x2) << " " << val << endl;
        }
        plotfile << endl;
    }
*/
    return 0;
}

