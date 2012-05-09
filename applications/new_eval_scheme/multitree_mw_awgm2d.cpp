#include <iostream>
#include <lawa/lawa.h>

using namespace std;
using namespace lawa;

typedef long double T;

typedef flens::GeMatrix<flens::FullStorage<T, cxxblas::ColMajor> >   DenseMatrixT;
typedef flens::DenseVector<flens::Array<T> >                         DenseVectorT;

typedef Basis<T,Orthogonal,Interval,Multi>                            PrimalBasis;
typedef PrimalBasis::RefinementBasis                                RefinementBasis;
typedef TensorBasis2D<Adaptive,PrimalBasis,PrimalBasis>             Basis2D;

typedef OptimizedH1Preconditioner2D<T,Basis2D>                      Preconditioner;

///  Underlying bilinear form
typedef RefinementBasis::LaplaceOperator1D                          RefinementLaplaceOp1D;

///  Local operator in 1d
typedef LocalOperator1D<PrimalBasis,PrimalBasis,
                        RefinementLaplaceOp1D>                      LocalOp1D;

typedef UniDirectionalLocalOperator<Index2D,XOne,LocalOp1D,
                                            NotXOne,Index1D>        UniDirectionalLocalOpXOne2D;
typedef UniDirectionalLocalOperator<Index2D,XTwo,LocalOp1D,
                                            NotXTwo,Index1D>        UniDirectionalLocalOpXTwo2D;

typedef CompoundLocalOperator<Index2D, UniDirectionalLocalOpXOne2D,
                              UniDirectionalLocalOpXTwo2D>          CompoundLocalOperator2D;

//Righthandsides definitions (separable)
typedef SeparableRHS2D<T,Basis2D >                                  SeparableRhsIntegral2D;

typedef SumOfTwoRHSIntegrals<T,Index2D,SeparableRhsIntegral2D,
                             SeparableRhsIntegral2D>                SumOfSeparableRhsIntegral2D;

typedef MultiTreeAWGM<Index2D,Basis2D,CompoundLocalOperator2D,
                      SumOfSeparableRhsIntegral2D,Preconditioner>   MultiTreeAWGM2D;

typedef IndexSet<Index1D>::const_iterator                           const_set1d_it;
typedef IndexSet<Index2D>::const_iterator                           const_set2d_it;
typedef Coefficients<Lexicographical,T,Index2D>::iterator           coeff2d_it;
typedef Coefficients<Lexicographical,T,Index2D>::const_iterator     const_coeff2d_it;

T u1(T x)   {    return x*x*(1-x)*(1-x); }

T u2(T y)   {    return y*y*(1-y)*(1-y); }

T du1(T x)  {    return 2*x*(1-x)*(1-x)-2*x*x*(1-x); }

T du2(T y)  {    return 2*y*(1-y)*(1-y)-2*y*y*(1-y); }

T ddu1(T x) {    return 2*(1-x)*(1-x) - 8*x*(1-x) + 2*x*x; }

T ddu2(T y) {    return 2*(1-y)*(1-y) - 8*y*(1-y) + 2*y*y; }

long double EnergyErrorSquared = 2.*(1.L/630.L * 2.L/105.L);

T f1(T x)   {   return -ddu1(x); }

T f2(T y)   {   return -ddu2(y); }

T sol(T x, T y) {   return u1(x) * u2(y); }

int main (int argc, char *argv[]) {

    cout.precision(20);
    if (argc!=4) {
        cout << "Usage: " << argv[0] << " d j0 J" << endl;
        return 0;
    }
    int d   = atoi(argv[1]);
    int j0  = atoi(argv[2]);
    int J  = atoi(argv[3]);
    T alpha = 0.7;
    T gamma = 0.005;
    T eps   = 1e-5;
    Timer time;

    /// Basis initialization
    //PrimalBasis       basis(d,d_,j0);
    PrimalBasis       basis(d,j0);
    basis.enforceBoundaryCondition<DirichletBC>();
    RefinementBasis  &refinementbasis = basis.refinementbasis;
    Basis2D basis2d(basis,basis);


    /// Operator initialization
    LocalOp1D                    localOp1D(basis,basis,refinementbasis.LaplaceOp1D);
    UniDirectionalLocalOpXOne2D  uniDirectionalOpXOne2D(localOp1D);
    UniDirectionalLocalOpXTwo2D  uniDirectionalOpXTwo2D(localOp1D);
    CompoundLocalOperator2D      localOp2D(uniDirectionalOpXOne2D,uniDirectionalOpXTwo2D);

    /// Initialization of preconditioner
    Preconditioner  Prec(basis2d,1.,1.,0.);

    /// Initialization of rhs
    DenseVectorT sing_pts_x, sing_pts_y;
    DenseMatrixT no_deltas, deltas_x, deltas_y;
    SeparableFunction2D<T> SepFunc1(f1, sing_pts_x, u2, sing_pts_y);
    SeparableFunction2D<T> SepFunc2(u1, sing_pts_x, f2, sing_pts_y);
    int order = 4+2*d;
    SeparableRhsIntegral2D      rhsintegral1(basis2d, SepFunc1, deltas_x, no_deltas, order);
    SeparableRhsIntegral2D      rhsintegral2(basis2d, SepFunc2, no_deltas, deltas_y, order);
    SumOfSeparableRhsIntegral2D F(rhsintegral1,rhsintegral2);

    /// Initialization of multi tree based adaptive wavelet Galerkin method
    MultiTreeAWGM2D multiTreeAWGM2D(basis2d, localOp2D, F, Prec);
    multiTreeAWGM2D.setParameters(alpha, gamma);

    Coefficients<Lexicographical,T,Index2D> u(SIZEHASHINDEX2D);
    getSparseGridVector(basis2d,u,0,(T)0.2);

    multiTreeAWGM2D.cg_solve(u, eps, "convfile.txt", 100, EnergyErrorSquared);

    plot2D<T,Basis2D,Preconditioner>(basis2d, u, Prec, sol, 0., 1., 0., 1., 0.1, "multiTreeAWGM_sol.txt");

    return 0;
}

