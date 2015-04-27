#include <iostream>
#include <lawa/lawa.h>
#include <applications/unbounded_domains/parameters/parameters.h>

using namespace std;
using namespace lawa;

typedef double T;

/// Thresh bound
const T thresh = 1e-15;

typedef flens::DenseVector<flens::Array<T> >                                DenseVectorT;

/// Basis definitions
typedef Basis<T,Orthogonal,Interval,Multi>                                  Basis1D;

/// Operator definitions
typedef AdaptiveHelmholtzOperatorOptimized1D<T,Orthogonal,Interval,Multi>   HelmholtzOp1D;


/// Iterators
typedef IndexSet<Index1D>::const_iterator                                   const_set_it;
typedef Coefficients<Lexicographical,T,Index1D>::const_iterator             const_coeff_it;

T
u1(T x) {
    return -exp(-100*(x-0.45)*(x-0.45));
    //return -x*x*(1-x)*(1-x);
}


T
f1(T x) {
    return (2*100-4*100*100*(x-0.45)*(x-0.45))*u1(x);
    //return 2*(1-x)*(1-x) - 8*x*(1-x) + 2*x*x;
}

int main (int argc, char *argv[]) {
    cout.precision(16);
    if (argc != 4) {
        cout << "usage " << argv[0] << " d j0 J" << endl;
        exit(1);
    }
    int d          =atoi(argv[1]);
    int j0         =atoi(argv[2]);   //minimal level for basis
    int J          =atoi(argv[3]);   //For each column index, add J row levels by lamdabTilde routine

    DenseVectorT singPts;
    Function<T> u_fct(u1,singPts);
    Function<T> f_fct(f1,singPts);

    Basis1D basis(d,j0);

    IntegralF<Gauss,Basis1D> integral_u(u_fct,basis);
    IntegralF<Gauss,Basis1D> integral_f(f_fct,basis);
    integral_u.quadrature.setOrder(20);
    integral_f.quadrature.setOrder(20);

    HelmholtzOp1D helmholtzOp1D(basis,0.);

    Coefficients<Lexicographical,T,Index1D> u, f;

    for (int k=basis.mra.rangeI(j0).firstIndex(); k<=basis.mra.rangeI(j0).lastIndex(); ++k) {
        Index1D index(j0,k,XBSpline);
        u[index] = 1./helmholtzOp1D.prec(index) * integral_u(j0,k,XBSpline,0);
        f[index] = helmholtzOp1D.prec(index) *    integral_f(j0,k,XBSpline,0);
    }
    for (int j=j0; j<=J; ++j) {
        for (int k=basis.rangeJ(j).firstIndex(); k<=basis.rangeJ(j).lastIndex(); ++k) {
            Index1D index(j,k,XWavelet);
            u[index] = 1./helmholtzOp1D.prec(index) * integral_u(j,k,XWavelet,0);
            f[index] = helmholtzOp1D.prec(index) *    integral_f(j,k,XWavelet,0);
        }
    }
    cout << "#supp u = " << u.size() << endl;
    cout << "#supp f = " << f.size() << endl;

    Coefficients<Lexicographical,T,Index1D> Au;
    ofstream file("conv_mw_apply1d_level.txt");
    for (int i=0; i<=J-j0-2; ++i) {
        Au = helmholtzOp1D.apply(u,i,J);
        int length1 = Au.size();
        Au -= f;
        T error =  Au.norm(2.);
        Coefficients<Lexicographical,T,Index1D> Au2;
        helmholtzOp1D.apply(u,error,Au2);
        int length2 = Au2.size();
        Au2 -= f;
        T error2 =  Au2.norm(2.);
        cout << i << " " << error << ", " << length1 << " || " << error2 << ", " << length2 << endl;
        file << i << " " << error << " " << length1 << " " << error2 << " " << length2 << endl;
    }
    file.close();

    return 0;
}

