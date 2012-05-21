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
typedef RHSWithPeaks1D<T,PrimalBasis>                               Rhs1D;
typedef AdaptiveSeparableRhs<T,Index2D,Rhs1D,Rhs1D >                AdaptiveSeparableRhsIntegral2D;
typedef CompoundRhs<T,Index2D,AdaptiveSeparableRhsIntegral2D,
                    AdaptiveSeparableRhsIntegral2D>                 CompoundRhsIntegral2D;

typedef MultiTreeAWGM<Index2D,Basis2D,CompoundLocalOperator2D,
                      CompoundRhsIntegral2D,Preconditioner>         MultiTreeAWGM2D;

typedef IndexSet<Index1D>::const_iterator                           const_set1d_it;
typedef IndexSet<Index2D>::const_iterator                           const_set2d_it;
typedef Coefficients<Lexicographical,T,Index2D>::iterator           coeff2d_it;
typedef Coefficients<Lexicographical,T,Index2D>::const_iterator     const_coeff2d_it;

int example = 2;

T u1(T x)   {    return 1.; }

T u2(T y)   {    return 1.; }

T du1(T x)  {    return 0.; }

T du2(T y)  {    return 0.; }

T ddu1(T x) {    return -10.;   }

T ddu2(T y) {    return -10.;   }

long double EnergyErrorSquared = 0.L;
/*
int example = 3;
T u1(T x)   {    return x*x*(1-x)*(1-x); }

T u2(T y)   {    return y*y*(1-y)*(1-y); }

T du1(T x)  {    return 2*x*(1-x)*(1-x)-2*x*x*(1-x); }

T du2(T y)  {    return 2*y*(1-y)*(1-y)-2*y*y*(1-y); }

T ddu1(T x) {    return 2*(1-x)*(1-x) - 8*x*(1-x) + 2*x*x; }

T ddu2(T y) {    return 2*(1-y)*(1-y) - 8*y*(1-y) + 2*y*y; }

long double EnergyErrorSquared = 2.*(1.L/630.L * 2.L/105.L);
*/
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
    int order = 4+2*d;

    Function<T>                    fct_u1(u1,sing_pts_x), fct_f1(f1,sing_pts_x);
    Function<T>                    fct_u2(u2,sing_pts_y), fct_f2(f2,sing_pts_y);
    RHSWithPeaks1D<T,PrimalBasis>  rhs_u1(basis, fct_u1, no_deltas, order);
    RHSWithPeaks1D<T,PrimalBasis>  rhs_f1(basis, fct_f1, no_deltas, order);
    RHSWithPeaks1D<T,PrimalBasis>  rhs_u2(basis, fct_u2, no_deltas, order);
    RHSWithPeaks1D<T,PrimalBasis>  rhs_f2(basis, fct_f2, no_deltas, order);
    Coefficients<Lexicographical,T,Index1D> rhs_u1_data(SIZEHASHINDEX1D),
                                            rhs_f1_data(SIZEHASHINDEX1D),
                                            rhs_u2_data(SIZEHASHINDEX1D),
                                            rhs_f2_data(SIZEHASHINDEX1D);
    AdaptiveSeparableRhsIntegral2D rhs1(rhs_f1, rhs_f1_data, rhs_u2, rhs_u2_data);
    AdaptiveSeparableRhsIntegral2D rhs2(rhs_u1, rhs_u1_data, rhs_f2, rhs_f2_data);
    CompoundRhsIntegral2D          F(rhs1,rhs2);

    /// Initialization of multi tree based adaptive wavelet Galerkin method
    MultiTreeAWGM2D multiTreeAWGM2D(basis2d, localOp2D, F, Prec);
    multiTreeAWGM2D.setParameters(alpha, gamma);


    Coefficients<Lexicographical,T,Index2D> u(SIZEHASHINDEX2D);
    getSparseGridVector(basis2d,u,0,(T)0.2);

    multiTreeAWGM2D.cg_solve(u, eps, "convfile.txt", 100, EnergyErrorSquared);

    plot2D<T,Basis2D,Preconditioner>(basis2d, u, Prec, sol, 0., 1., 0., 1., 0.1, "multiTreeAWGM_sol");

    return 0;
}

