#include <iostream>
#include <lawa/lawa.h>

using namespace std;
using namespace lawa;

typedef long double T;

typedef flens::GeMatrix<flens::FullStorage<T, cxxblas::ColMajor> >  DenseMatrixT;
typedef flens::DenseVector<flens::Array<T> >                        DenseVectorT;

typedef Basis<T,Orthogonal,Interval,Multi>                          PrimalBasis;
typedef PrimalBasis::RefinementBasis                                RefinementBasis;
typedef TensorBasis2D<Adaptive,PrimalBasis,PrimalBasis>             Basis2D;

typedef OptimizedH1Preconditioner2D<T,Basis2D>                      Preconditioner;

///  Underlying bilinear form
typedef RefinementBasis::LaplaceOperator1D                          RefinementLaplaceOp1D;
typedef AdaptiveLaplaceOperator1D<T,Orthogonal,Interval,Multi>      LaplaceOp1D;


///  Local operator in 1d
typedef LocalOperator1D<PrimalBasis,PrimalBasis,
                        RefinementLaplaceOp1D,LaplaceOp1D>          LocalOp1D;

typedef UniDirectionalLocalOperator<Index2D,XOne,LocalOp1D,
                                            NotXOne,Index1D>        UniDirectionalLocalOpXOne2D;
typedef UniDirectionalLocalOperator<Index2D,XTwo,LocalOp1D,
                                            NotXTwo,Index1D>        UniDirectionalLocalOpXTwo2D;

typedef CompoundLocalOperator<Index2D, UniDirectionalLocalOpXOne2D,
                              UniDirectionalLocalOpXTwo2D>          CompoundLocalOperator2D;

typedef ThetaTimeStepLocalOperator<Index2D, CompoundLocalOperator2D> ThetaTimeStepLocalOperator2D;

//Righthandsides definitions (separable)
typedef RHSWithPeaks1D<T,PrimalBasis>                               Rhs1D;
typedef AdaptiveSeparableRhs<T,Index2D,Rhs1D,Rhs1D >                AdaptiveSeparableRhsIntegral2D;
typedef ThetaTimeStepSeparableRHS<T,Index2D,
                                  AdaptiveSeparableRhsIntegral2D,
                                  ThetaTimeStepLocalOperator2D>   ThetaTimeStepRhs2d;

typedef IndexSet<Index1D>::const_iterator                           const_set1d_it;
typedef IndexSet<Index2D>::const_iterator                           const_set2d_it;
typedef Coefficients<Lexicographical,T,Index1D>::iterator           coeff1d_it;
typedef Coefficients<Lexicographical,T,Index2D>::iterator           coeff2d_it;
typedef Coefficients<Lexicographical,T,Index1D>::const_iterator     const_coeff1d_it;
typedef Coefficients<Lexicographical,T,Index2D>::const_iterator     const_coeff2d_it;

T time_point = 0.;
T theta = 0.5;
T timestep = 0.25;
size_t hms = 49157;

T sol(T t, T x, T y)   {    return (std::exp(t)-1.) * std::sin(M_PI*x) * std::sin(2*M_PI*y) ; }

T f_t(T t) {  return  (5.*M_PI*M_PI+1)*std::exp(t) - 5.*M_PI*M_PI; }
T f1(T x) {  return  std::sin(M_PI*x); }
T f2(T y) {  return  std::sin(2.*M_PI*y); }


void
multiply_by_Prec(const Preconditioner& P, Coefficients<Lexicographical,T,Index2D> &v);

void
multiply_by_invPrec(const Preconditioner& P, Coefficients<Lexicographical,T,Index2D> &v);

void
mycg(ThetaTimeStepLocalOperator2D& A, const Preconditioner& P,
     const Coefficients<Lexicographical,T,Index2D>& f, Coefficients<Lexicographical,T,Index2D> &u,
     T tol);

T
plot_current_solution(T time, const Basis2D& basis2d, const Coefficients<Lexicographical,T,Index2D> &u);

int main (int argc, char *argv[]) {

    cout.precision(20);
    if (argc!=4) {
        cout << "Usage: " << argv[0] << " d j0 J" << endl;
        return 0;
    }

    int d   = atoi(argv[1]);
    int j0  = atoi(argv[2]);
    int j  = atoi(argv[3]);
    T alpha = 0.7;
    T gamma = 0.1;
    const char* residualType = "standard";
    const char* treeType = "sparsetree"; //"gradedtree";
    bool IsMW = true;

    T eps   = 1e-5;
    Timer time;


    /// Basis initialization
    //PrimalBasis       basis(d,d_,j0);
    PrimalBasis       basis(d,j0);
    basis.enforceBoundaryCondition<DirichletBC>();
    RefinementBasis  &refinementbasis = basis.refinementbasis;
    Basis2D basis2d(basis,basis);

    /// Operator initialization
    LaplaceOp1D                  laplaceOp1D(basis);
    LocalOp1D                    localOp1D(basis,basis,refinementbasis.LaplaceOp1D,laplaceOp1D);
    UniDirectionalLocalOpXOne2D  uniDirectionalOpXOne2D(localOp1D);
    UniDirectionalLocalOpXTwo2D  uniDirectionalOpXTwo2D(localOp1D);
    CompoundLocalOperator2D      localOp2D(uniDirectionalOpXOne2D,uniDirectionalOpXTwo2D);
    ThetaTimeStepLocalOperator2D localThetaTimeStepOp2D(theta,timestep,localOp2D);

    /// Initialization of preconditioner
    Preconditioner  Prec(basis2d,theta*timestep, theta*timestep, 1.);

    /// Initialization of rhs
    DenseVectorT sing_pts_t, sing_pts_x, sing_pts_y;
    DenseMatrixT no_deltas, deltas_x, deltas_y;
    int order = 7;
    Function<T>                    fct_f_t(f_t,sing_pts_t);
    Function<T>                    fct_f1(f1,sing_pts_x), fct_f2(f2,sing_pts_y);
    RHSWithPeaks1D<T,PrimalBasis>  rhs_f1(basis, fct_f1, no_deltas, order);
    RHSWithPeaks1D<T,PrimalBasis>  rhs_f2(basis, fct_f2, no_deltas, order);
    Coefficients<Lexicographical,T,Index1D> rhs_f1_data(SIZEHASHINDEX1D),
                                            rhs_f2_data(SIZEHASHINDEX1D);
    AdaptiveSeparableRhsIntegral2D rhs(rhs_f1, rhs_f1_data, rhs_f2, rhs_f2_data);
    ThetaTimeStepRhs2d thetatimestep_rhs(fct_f_t,rhs,localThetaTimeStepOp2D);

    stringstream convfilename;
    convfilename << "conv_theta_" << theta << "_d_" << d << "_j_" << j << ".dat";
    ofstream convfile(convfilename.str().c_str());

    for (int timesteps=1; timesteps<=32; timesteps*=2) {
        timestep = 1./timesteps;
        std::cerr << "timestep = " << timestep << std::endl;
        T max_error = 0.;
        Coefficients<Lexicographical,T,Index2D> u_k(hms), u_kP1(hms);
        Coefficients<Lexicographical,T,Index2D> propagated_u_k(hms);
        Coefficients<Lexicographical,T,Index2D> f(hms);

        getSparseGridVector(basis2d, propagated_u_k, j, 0.L);
        getSparseGridVector(basis2d, u_k, j, 0.L);
        getSparseGridVector(basis2d, u_kP1, j, 0.L);
        getSparseGridVector(basis2d, f, j, 0.L);

        for (int i=0; i<timesteps; ++i) {
            Prec.setParameters(theta*timestep, theta*timestep, 1.);
            localThetaTimeStepOp2D.setThetaTimeStepParameters(theta,timestep);
            T time = (i+1)*timestep;
            std::cerr << "Current discrete time point: " << time << std::endl;
            f.setToZero();
            propagated_u_k.setToZero();
            if (theta!=1) {
                localThetaTimeStepOp2D.evalA(u_k, propagated_u_k, "galerkin");
                propagated_u_k *= (theta-1.)*timestep;
            }
            propagated_u_k += u_k;
            thetatimestep_rhs.setThetaTimeStepParameters(theta, timestep, time, propagated_u_k);
            for (coeff2d_it it=f.begin(); it!=f.end(); ++it) {
                (*it).second = thetatimestep_rhs((*it).first);
            }


            multiply_by_Prec(Prec, f);
            mycg(localThetaTimeStepOp2D, Prec, f, u_kP1, 1e-6);
            u_k = u_kP1;
            multiply_by_Prec(Prec, u_k);
            cout << "u = " << u_k << endl;
        }
        max_error = plot_current_solution(0.25, basis2d, u_k);
        convfile << timestep << " " << max_error << endl;
        //cout << "Hit enter" << endl;
        //getchar();
    }

    return 0;
}


void
multiply_by_Prec(const Preconditioner& P, Coefficients<Lexicographical,T,Index2D> &v)
{
    for (coeff2d_it it=v.begin(); it!=v.end(); ++it) {
        (*it).second *= P[(*it).first];
    }
}

void
multiply_by_invPrec(const Preconditioner& P, Coefficients<Lexicographical,T,Index2D> &v)
{
    for (coeff2d_it it=v.begin(); it!=v.end(); ++it) {
        (*it).second *= 1./P[(*it).first];
    }
}

void
mycg(ThetaTimeStepLocalOperator2D& A, const Preconditioner& P,
   const Coefficients<Lexicographical,T,Index2D>& f, Coefficients<Lexicographical,T,Index2D> &u,
   T tol)
{
    Coefficients<Lexicographical,T,Index2D> p(hms), r(hms), Ap(hms);
    p = u; r = u; Ap = u;
    p.setToZero();
    r.setToZero();
    Ap.setToZero();

    A.eval(u,r,P,"galerkin");
    r -= f;
    p = r;
    p *= (T)(-1.);
    T cg_rNormSquare = r*r;
    int maxIterations=100;
    int cg_iter=0;
    for (cg_iter=0; cg_iter<maxIterations; ++cg_iter) {
        if (std::sqrt(cg_rNormSquare)<=tol) {
            std::cerr << "      CG stopped after " << cg_iter << " iterations with error "
                      << sqrt(cg_rNormSquare) << std::endl;
            break;
        }
        //std::cerr << "    Iteration " << cg_iter+1 << std::endl;
        Ap.setToZero();
        A.eval(p,Ap,P,"galerkin");

        T pAp = p * Ap;
        T alpha = cg_rNormSquare/pAp;
        p *= alpha;
        u += p;
        p *= (T)1./alpha;
        Ap *= alpha;
        r += Ap;

        T cg_rNormSquarePrev = cg_rNormSquare;
        cg_rNormSquare = r*r;
        //std::cerr << "      Current error in cg: " << std::sqrt(cg_rNormSquare) << std::endl;
        T beta = cg_rNormSquare/cg_rNormSquarePrev;
        p *= beta;
        p -= r;
    }
}

T
plot_current_solution(T time, const Basis2D& basis2d, const Coefficients<Lexicographical,T,Index2D> &u)
{
    stringstream plotfilename;
    plotfilename << "theta_sol_" << time << ".dat";
    ofstream plotfile(plotfilename.str().c_str());

    T max_error = 0.;
    for (T x=0.; x<=1.; x+=0.01) {
        for (T y=0.; y<=1.; y+=0.01) {
            T val = 0.;
            for (const_coeff2d_it it=u.begin(); it!=u.end(); ++it) {
                XType xtype_x = (*it).first.index1.xtype, xtype_y = (*it).first.index2.xtype;
                int   j_x = (*it).first.index1.j, j_y = (*it).first.index2.j;
                long  k_x = (*it).first.index1.k, k_y = (*it).first.index2.k;
                val += (*it).second * basis2d.first.generator(xtype_x)(x,j_x,k_x,0) * basis2d.second.generator(xtype_y)(y,j_y,k_y,0);
                //std::cerr << "val = " << val << std::endl;
            }
            plotfile << x << " " << y << " " << sol(time,x,y) << " " << val << endl;
            max_error = std::max(max_error, fabs(val-sol(time,x,y)));
        }
        plotfile << endl;
    }
    std::cerr << "Maximum error: " << max_error << std::endl;
    return max_error;
}

