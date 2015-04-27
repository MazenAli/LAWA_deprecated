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
typedef TensorBasis2D<Adaptive,PrimalBasis,PrimalBasis>                     Basis2D;

/// Operator definitions
typedef AdaptiveHelmholtzOperatorOptimized2D<T,Orthogonal,Interval,Multi,
                                               Orthogonal,Interval,Multi>   HelmholtzOp2D;

typedef OptimizedH1Preconditioner2D<T,Basis2D>                              Preconditioner;

///  Underlying bilinear form
typedef RefinementBasis::LaplaceOperator1D                                  RefinementLaplaceOp1D;
typedef AdaptiveLaplaceOperator1D<T,Orthogonal,Interval,Multi>              LaplaceOp1D;

///  Local operator in 1d
typedef LocalOperator1D<PrimalBasis,PrimalBasis,RefinementLaplaceOp1D,
                        LaplaceOp1D>                                        LocalOp1D;

typedef UniDirectionalLocalOperator<Index2D,XOne,LocalOp1D,
                                            NotXOne,Index1D>                UniDirectionalLocalOpXOne2D;
typedef UniDirectionalLocalOperator<Index2D,XTwo,LocalOp1D,
                                            NotXTwo,Index1D>                UniDirectionalLocalOpXTwo2D;

typedef CompoundLocalOperator<Index2D, UniDirectionalLocalOpXOne2D,
                              UniDirectionalLocalOpXTwo2D>          CompoundLocalOperator2D;

/// Iterators
typedef Coefficients<Lexicographical,T,Index1D>::iterator                   coeff1d_it;
typedef Coefficients<Lexicographical,T,Index2D>::iterator                   coeff2d_it;

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
    Basis2D basis2d(basis,basis);

    HelmholtzOp2D helmholtzOp2D(basis2d,0.);

    /// Operator initialization
    LaplaceOp1D                  laplaceOp1D(basis);
    LocalOp1D                    localOp1D(basis,basis,refinementbasis.LaplaceOp1D,laplaceOp1D);
    UniDirectionalLocalOpXOne2D  uniDirectionalOpXOne2D(localOp1D);
    UniDirectionalLocalOpXTwo2D  uniDirectionalOpXTwo2D(localOp1D);
    CompoundLocalOperator2D      localOp2D(uniDirectionalOpXOne2D,uniDirectionalOpXTwo2D);

    /// Initialization of preconditioner
    Preconditioner  Prec(basis2d,1.,1.,0.);

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

    Coefficients<Lexicographical,T,Index2D> u, f;
    getSparseGridVector(basis2d, u, J, 0.);
    getSparseGridVector(basis2d, f, J, 0.);

    for (coeff2d_it it=u.begin(); it!=u.end(); ++it) {
        (*it).second = help_u[(*it).first.index1] * help_u[(*it).first.index2]
                     * 1./helmholtzOp2D.prec((*it).first);
    }
    for (coeff2d_it it=f.begin(); it!=f.end(); ++it) {
        (*it).second  = help_u[(*it).first.index1] * help_f[(*it).first.index2];
        (*it).second += help_f[(*it).first.index1] * help_u[(*it).first.index2];
        (*it).second *= helmholtzOp2D.prec((*it).first);
    }
    cout << "#supp u = " << u.size() << endl;
    cout << "#supp f = " << f.size() << endl;

    Coefficients<Lexicographical,T,Index2D> Au;
    ofstream file("conv_mw_apply2d_level.txt");
    Timer time;
    for (int i=0; i<=J; ++i) {
        time.start();
        Au = helmholtzOp2D.apply(u,i,J);
        time.stop();
        T time_apply1 = time.elapsed();
        int length1 = Au.size();
        Au -= f;
        T error =  Au.norm(2.);

        Au.clear();
        time.start();
        helmholtzOp2D.apply(u,error,Au);
        time.stop();
        T time_apply2 = time.elapsed();
        int length2 = Au.size();
        Au -= f;
        T error2 =  Au.norm(2.);

        Au.clear();
        time.start();
        localOp2D.apply(u, Au, Prec, error);
        time.stop();
        T time_apply3 = time.elapsed();
        int length3 = Au.size();
        Au -= f;
        T error3 =  Au.norm(2.);



        cout << i << " " << error << " " << length1 << " " << time_apply1 << " || "
                         << error2 << " " << length2 << " " << time_apply2 << " || "
                         << error3 << " " << length3 << " " << time_apply3 << endl;
        file << i << " " << error << " " << length1 << " " << error2 << " " << length2 << endl;
    }
    file.close();
    return 0;

}

