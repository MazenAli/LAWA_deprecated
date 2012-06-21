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
typedef Basis<T,Orthogonal,Interval,Multi>                                  PrimalBasis;
typedef PrimalBasis::RefinementBasis                                        RefinementBasis;
typedef TensorBasis3D<Adaptive,PrimalBasis,PrimalBasis,PrimalBasis>         Basis3D;

/// Operator definitions
typedef OptimizedH1Preconditioner3D<T,Basis3D>                              Preconditioner;

///  Underlying bilinear form
typedef RefinementBasis::LaplaceOperator1D                                  RefinementLaplaceOp1D;
typedef AdaptiveLaplaceOperator1D<T,Orthogonal,Interval,Multi>              LaplaceOp1D;

///  Local operator in 1d
typedef LocalOperator1D<PrimalBasis,PrimalBasis,RefinementLaplaceOp1D,
                        LaplaceOp1D>                                        LocalOp1D;

typedef UniDirectionalLocalOperator<Index3D,XOne,LocalOp1D,
                                            NotXOne,Index2D>                UniDirectionalLocalOpXOne3D;
typedef UniDirectionalLocalOperator<Index3D,XTwo,LocalOp1D,
                                            NotXTwo,Index2D>                UniDirectionalLocalOpXTwo3D;
typedef UniDirectionalLocalOperator<Index3D,XThree,LocalOp1D,
                                            NotXThree,Index2D>              UniDirectionalLocalOpXThree3D;

typedef CompoundLocalOperator<Index3D, UniDirectionalLocalOpXOne3D,
                                       UniDirectionalLocalOpXTwo3D,
                                       UniDirectionalLocalOpXThree3D>       CompoundLocalOperator3D;

/// Iterators
typedef Coefficients<Lexicographical,T,Index1D>::iterator                   coeff1d_it;
typedef Coefficients<Lexicographical,T,Index2D>::iterator                   coeff2d_it;
typedef Coefficients<Lexicographical,T,Index3D>::iterator                   coeff3d_it;

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

    PrimalBasis basis(d,j0);
    RefinementBasis  &refinementbasis = basis.refinementbasis;
    Basis3D basis3d(basis,basis,basis);

    /// Operator initialization
    LaplaceOp1D                     laplaceOp1D(basis);
    LocalOp1D                       localOp1D(basis,basis,refinementbasis.LaplaceOp1D,laplaceOp1D);
    UniDirectionalLocalOpXOne3D     uniDirectionalOpXOne3D(localOp1D);
    UniDirectionalLocalOpXTwo3D     uniDirectionalOpXTwo3D(localOp1D);
    UniDirectionalLocalOpXThree3D   uniDirectionalOpXThree3D(localOp1D);
    CompoundLocalOperator3D         localOp3D(uniDirectionalOpXOne3D,uniDirectionalOpXTwo3D,
                                              uniDirectionalOpXThree3D);

    /// Initialization of preconditioner
    Preconditioner  Prec(basis3d,1.,1.,1.,0.);

    Coefficients<Lexicographical,T,Index1D> help_u, help_f;

    IntegralF<Gauss,PrimalBasis> integral_u(u_fct,basis);
    IntegralF<Gauss,PrimalBasis> integral_f(f_fct,basis);
    integral_u.quadrature.setOrder(20);
    integral_f.quadrature.setOrder(20);
    for (int k=basis.mra.rangeI(j0).firstIndex(); k<=basis.mra.rangeI(j0).lastIndex(); ++k) {
        Index1D index(j0,k,XBSpline);
        help_u[index] = integral_u(j0,k,XBSpline,0);
        help_f[index] = integral_f(j0,k,XBSpline,0);
    }
    for (int j=j0; j<=J; ++j) {
        for (int k=basis.rangeJ(j).firstIndex(); k<=basis.rangeJ(j).lastIndex(); ++k) {
            Index1D index(j,k,XWavelet);
            help_u[index] = integral_u(j,k,XWavelet,0);
            help_f[index] = integral_f(j,k,XWavelet,0);
        }
    }
    cout << "#supp help_u = " << help_u.size() << endl;
    cout << "#supp help_f = " << help_f.size() << endl;

    Coefficients<Lexicographical,T,Index3D> u(SIZELARGEHASHINDEX2D), f(SIZELARGEHASHINDEX2D);
    getSparseGridVector(basis3d, u, J, 0.);
    getSparseGridVector(basis3d, f, J, 0.);


    for (coeff3d_it it=u.begin(); it!=u.end(); ++it) {
        (*it).second = help_u[(*it).first.index1] * help_u[(*it).first.index2]
                     * help_u[(*it).first.index3] * 1./Prec((*it).first);
    }
    for (coeff3d_it it=f.begin(); it!=f.end(); ++it) {
        (*it).second  = help_f[(*it).first.index1] * help_u[(*it).first.index2] * help_u[(*it).first.index3];
        (*it).second += help_u[(*it).first.index1] * help_f[(*it).first.index2] * help_u[(*it).first.index3];
        (*it).second += help_u[(*it).first.index1] * help_u[(*it).first.index2] * help_f[(*it).first.index3];
        (*it).second *= Prec((*it).first);

    }
    cout << "#supp u = " << u.size() << endl;
    cout << "#supp f = " << f.size() << endl;

    T tol = f.norm();
    //cout << "u = " << u << endl;

    while (tol>1e-3) {
        Coefficients<Lexicographical,T,Index3D> Au(SIZELARGEHASHINDEX2D);
        localOp3D.apply(u, Au, Prec, tol);
        Au -= f;
        cout << "Target accuracy: tol = " << tol << ", attained accuracy: " << Au.norm(2.) << endl;
        tol *= 0.5;
    }

    return 0;
}
