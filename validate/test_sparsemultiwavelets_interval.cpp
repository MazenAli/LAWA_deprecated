#include <iostream>
#include <lawa/lawa.h>

using namespace std;
using namespace lawa;

typedef double T;

typedef flens::DenseVector<flens::Array<T> > DenseVectorT;

typedef Basis<T,Primal,Interval,SparseMulti> Basis1D;

typedef Integral<Gauss,Basis1D,Basis1D>     Integral1D;
typedef IntegralF<Gauss,Basis1D>            IntegralF1D;

T
p0(T x) {
    return 1.;
}

T
p1(T x) {
    return x;
}

T
p2(T x) {
    return x*x;
}

T
p3(T x) {
    return x*x*x;
}

Range<long>
getRange(const Basis1D& basis, int j, XType type);

int main (int argc, char *argv[]) {

    cout.precision(16);
    int d=4;
    int j0=1;
    int j=1;

    Basis1D basis(d,j0);
    basis.enforceBoundaryCondition<DirichletBC>();

    Integral1D integral(basis,basis);

    DenseVectorT singPts;
    Function<T> p0Fct(p0,singPts);
    Function<T> p1Fct(p1,singPts);
    Function<T> p2Fct(p2,singPts);
    Function<T> p3Fct(p3,singPts);

    IntegralF1D integralp0(p0Fct,basis);
    IntegralF1D integralp1(p1Fct,basis);
    IntegralF1D integralp2(p2Fct,basis);
    IntegralF1D integralp3(p3Fct,basis);

    Range<long> waveletrange = getRange(basis, j, XWavelet);
    for(long k=waveletrange.firstIndex(); k<=waveletrange.lastIndex(); ++k) {
        ofstream file("sparsemulti_wavelet.dat");
        cout << basis.psi.support(j,k) << endl;
        cout << basis.psi.singularSupport(j,k) << endl;
        cout << "<psi, psi> = " << integral(j,k,XWavelet,0,j,k,XWavelet,0) << endl;
        cout << "<p0, psi> = " << integralp0(j,k,XWavelet,0) << endl;
        cout << "<p1, psi> = " << integralp1(j,k,XWavelet,0) << endl;
        cout << "<p2, psi> = " << integralp2(j,k,XWavelet,0) << endl;
        cout << "<p3, psi> = " << integralp3(j,k,XWavelet,0) << endl;
        T l1=basis.psi.support(j,k).l1;
        T l2=basis.psi.support(j,k).l2;
        for (T x=l1; x<=l2; x+=pow2i<T>(-9)) {
            file << x << " " << basis.psi(x,j,k,0) << " " << basis.psi(x,j,k,1)  << endl;
        }
        file.close();
        getchar();
    }

    Range<long> scalingrange = getRange(basis, j, XBSpline);
    for (long k=scalingrange.firstIndex(); k<=scalingrange.lastIndex(); ++k) {
        ofstream file("sparsemulti_scaling.dat");

        cout << basis.mra.phi.support(j,k) << endl;
        cout << basis.mra.phi.singularSupport(j,k) << endl;
        cout << "<phi, phi> = " << integral(j,k,XBSpline,0,j,k,XBSpline,0) << endl;
        for (T x=-2.; x<=2.; x+=pow2i<T>(-8)) {
            file << x << " " << basis.mra.phi(x,j,k,0) << " " << basis.mra.phi(x,j,k,1) <<  endl;
        }
        file.close();
        getchar();
    }

    return 0;
}

Range<long>
getRange(const Basis1D& basis, int j, XType type)
{
    if (type==XBSpline) {
        return basis.mra.rangeI(j);
    }
    else {
        return basis.rangeJ(j);
    }
}
