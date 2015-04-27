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

typedef MultiTreeAWGM<Index2D,Basis2D,ThetaTimeStepLocalOperator2D,
                      ThetaTimeStepRhs2d,Preconditioner>            ThetaTimeStepMultiTreeAWGM2D;

typedef MultiTreeAWGM<Index2D,Basis2D,
                      ThetaTimeStepLocalOperator2D,
                      AdaptiveSeparableRhsIntegral2D,
                      NoPreconditioner<T,Index2D> >                 ApproxL2AWGM2D;

typedef ThetaSchemeAWGM<Index2D, ThetaTimeStepMultiTreeAWGM2D>      ThetaSchemeMultiTreeAWGM2D;

typedef IndexSet<Index1D>::const_iterator                           const_set1d_it;
typedef IndexSet<Index2D>::const_iterator                           const_set2d_it;
typedef Coefficients<Lexicographical,T,Index1D>::iterator           coeff1d_it;
typedef Coefficients<Lexicographical,T,Index2D>::iterator           coeff2d_it;
typedef Coefficients<Lexicographical,T,Index1D>::const_iterator     const_coeff1d_it;
typedef Coefficients<Lexicographical,T,Index2D>::const_iterator     const_coeff2d_it;


int example = 1;
T sol(T t, T x, T y)   {    return (std::exp(t)) * std::sin(M_PI*x) * std::sin(2*M_PI*y) ; }

T f_t(T t)  {  return  (5.*M_PI*M_PI+1)*std::exp(t); }
T f1(T x)   {  return  std::sin(M_PI*x); }
T f2(T y)   {  return  std::sin(2.*M_PI*y); }
T u1(T x)   {  return  std::sin(M_PI*x); }
T u2(T y)   {  return  std::sin(2.*M_PI*y); }

/*
int example = 2;
T sol(T t, T x, T y)   {    return std::pow(t,(T)0.75) * std::sin(M_PI*x) * std::sin(2*M_PI*y) ; }

T f_t(T t)  { if (t<=0) std::cerr << "called for t = " << t << std::endl;
    return  0.75*std::pow(t,(T)-0.25)+(5.*M_PI*M_PI)*std::pow(t,(T)0.75); }
T f1(T x)   {  return  std::sin(M_PI*x); }
T f2(T y)   {  return  std::sin(2.*M_PI*y); }
T u1(T x)   {  return  0.; }
T u2(T x)   {  return  0.; }
*/


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
    T gamma = 0.025;
    const char* residualType = "standard";
    const char* treeType = "sparsetree"; //"gradedtree";
    bool IsMW = true;

    T theta = 0.5;
    T timestep = 0.25;
    T timestep_eps = 1e-4;
    int maxiterations = 50;  T init_cgtol = 1e-2;
    //int maxiterations =  1;  T init_cgtol = 1e-8;   // use maxiterations = 1 for "pure" sparse grid computation
    //int numOfTimesteps = 4;

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
    NoPreconditioner<T, Index2D> NoPrec;
    Preconditioner  Prec(basis2d,theta*timestep, theta*timestep, 1.);

    /// Initialization of integrals for initial condition and rhs
    DenseVectorT sing_pts_t, sing_pts_x, sing_pts_y;
    DenseMatrixT no_deltas, deltas_x, deltas_y;
    int order = 7;
    Function<T>                    fct_f_t(f_t,sing_pts_t);
    Function<T>                    fct_f1(f1,sing_pts_x), fct_f2(f2,sing_pts_y);
    Function<T>                    fct_u1(u1,sing_pts_x), fct_u2(u2,sing_pts_y);
    RHSWithPeaks1D<T,PrimalBasis>  rhs_f1(basis, fct_f1, no_deltas, order);
    RHSWithPeaks1D<T,PrimalBasis>  rhs_f2(basis, fct_f2, no_deltas, order);
    RHSWithPeaks1D<T,PrimalBasis>  rhs_u1(basis, fct_u1, no_deltas, order);
    RHSWithPeaks1D<T,PrimalBasis>  rhs_u2(basis, fct_u2, no_deltas, order);
    Coefficients<Lexicographical,T,Index1D> rhs_f1_data(SIZEHASHINDEX1D),
                                            rhs_f2_data(SIZEHASHINDEX1D),
                                            rhs_u1_data(SIZEHASHINDEX1D),
                                            rhs_u2_data(SIZEHASHINDEX1D);
    AdaptiveSeparableRhsIntegral2D F_rhs(rhs_f1, rhs_f1_data, rhs_f2, rhs_f2_data);
    AdaptiveSeparableRhsIntegral2D F_u0(rhs_u1, rhs_u1_data, rhs_u2, rhs_u2_data);

    ThetaTimeStepRhs2d thetatimestep_F(fct_f_t,F_rhs,localThetaTimeStepOp2D);

    /// Initialization of multi tree based adaptive wavelet Galerkin method
    ApproxL2AWGM2D approxL2_solver(basis2d, localThetaTimeStepOp2D, F_u0, NoPrec);
    ThetaTimeStepMultiTreeAWGM2D thetatimestep_solver(basis2d, localThetaTimeStepOp2D,
                                                      thetatimestep_F, Prec);
    approxL2_solver.setParameters(alpha, gamma, residualType, treeType, IsMW, false);
    thetatimestep_solver.setParameters(alpha, gamma, residualType, treeType, IsMW, false);

    ThetaSchemeMultiTreeAWGM2D thetascheme(thetatimestep_solver);

    stringstream convfilename;
    if (maxiterations == 1) {
        convfilename << "conv_ex_" << example << "_theta_" << theta << "_d_" << d << "_j_" << j << ".dat";
    }
    else {
        convfilename << "_conv_ex_" << example << "_theta_" << theta << "_d_" << d << "_eps_" << timestep_eps << ".dat";
    }
    ofstream convfile(convfilename.str().c_str());

    for (int numOfTimesteps=1; numOfTimesteps<=64; numOfTimesteps*=2) {
        timestep = 1./numOfTimesteps;
        //thetascheme.setParameters(theta, timestep, numOfTimesteps, timestep_eps, maxiterations,
        //                          init_cgtol);

        size_t hms = 49157;
        Coefficients<Lexicographical,T,Index2D> u(hms);
        getSparseGridVector(basis2d, u, j, 0.L);

        //approxL2_solver.approxL2(u, timestep_eps, maxiterations);
        //if (example = 1) {
        //    for (coeff2d_it it=u.begin(); it!=u.end(); ++it) { (*it).second = rhs((*it).first); }
        //}

        //thetascheme.solve(u);
        cerr << "Warning: solve commented out, incompatable" << endl;
        T maxerror = plot_current_solution(numOfTimesteps*timestep, basis2d, u);
        convfile << timestep << " " << maxerror << endl;

    }
    convfile.close();

    return 0;
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

/*
 *         Coefficients<Lexicographical,T,Index2D> u_test(hms), dummy(hms);
        getSparseGridVector(basis2d, u_test, 0, 0.L);
        getSparseGridVector(basis2d, dummy, 0, 0.L);
        thetatimestep_solver.F.setThetaTimeStepParameters(1.,1.,0.,&dummy);
        thetatimestep_solver.approxL2(u, timestep_eps, 20);
        plot_current_solution(0., basis2d, u);
 *
 */
