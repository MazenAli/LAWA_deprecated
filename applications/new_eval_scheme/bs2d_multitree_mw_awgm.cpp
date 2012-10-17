#include <iostream>
#include <lawa/lawa.h>

#include <applications/finance/options/options.h>
#include <applications/finance/operators/operators.h>
#include <applications/finance/processes/processes.h>


using namespace std;
using namespace lawa;

typedef long double T;


const OptionType1D optiontype = Put;
T strike = 1.;
T maturity = 1.;
OptionParameters1D<T,Put> optionparameters(strike, maturity, false);

T r = 0.;
T sigma1 = 0.3, sigma2 = 0.1;
const ProcessType1D  processtype  = BlackScholes;
ProcessParameters1D<T,BlackScholes>   processparameters1(r, sigma1);
ProcessParameters1D<T,BlackScholes>   processparameters2(r, sigma2);

T weight1 = 1.;
T weight2 = 1.;

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
                                  AdaptiveSeparableRhsIntegral2D>   ThetaTimeStepRhs2d;
typedef CompoundRhs<T,Index2D,AdaptiveSeparableRhsIntegral2D,
                    AdaptiveSeparableRhsIntegral2D>                 CompoundRhsIntegral2D;

typedef MultiTreeAWGM<Index2D,Basis2D,ThetaTimeStepLocalOperator2D,
                      ThetaTimeStepRhs2d,Preconditioner>            ThetaTimeStepMultiTreeAWGM2D;

typedef MultiTreeAWGM<Index2D,Basis2D,
                      ThetaTimeStepLocalOperator2D,
                      CompoundRhsIntegral2D,
                      NoPreconditioner<T,Index2D> >                 ApproxL2AWGM2D;

typedef ThetaSchemeAWGM<Index2D, ThetaTimeStepMultiTreeAWGM2D>      ThetaSchemeMultiTreeAWGM2D;

typedef IndexSet<Index1D>::const_iterator                           const_set1d_it;
typedef IndexSet<Index2D>::const_iterator                           const_set2d_it;
typedef Coefficients<Lexicographical,T,Index1D>::iterator           coeff1d_it;
typedef Coefficients<Lexicographical,T,Index2D>::iterator           coeff2d_it;
typedef Coefficients<Lexicographical,T,Index1D>::const_iterator     const_coeff1d_it;
typedef Coefficients<Lexicographical,T,Index2D>::const_iterator     const_coeff2d_it;

T a = -4.;
T b =  4.;

T decaywidth = 1.;

T alpha = 1.L  / 100.L;
T beta  = 99.L / 100.L;

T trunc(T x) {
    if (x <= alpha)                 return 1.;
    else if (alpha < x && x < beta) {
        T numerator = std::exp(-1.L/((beta-x)*(beta-x)));
        T denominator = std::exp(-1.L/((x-alpha)*(x-alpha))) + numerator;
        return numerator / denominator;
    }
    else return 0.;
}

T sol(T t, T x, T y) { return 1.; }

T dom_map(T x) { return (b-a)*x + a; }

T u0(T x, T y)   {  return   (weight1*max(1.-exp(dom_map(x)),0.L) + weight2*max(1.-exp(dom_map(y)),0.L))
                           * trunc(a+decaywidth-dom_map(x)) * trunc(a+decaywidth-dom_map(y))
                           * trunc(decaywidth-b+dom_map(x)) * trunc(decaywidth-b+dom_map(y)); }
T u1_x(T x)      {  return  weight1*max(1.-exp(dom_map(x)),0.L) * trunc(a+decaywidth-dom_map(x)) * trunc(decaywidth-b+dom_map(x)); }
T u1_y(T y)      {  return  1. * trunc(a+decaywidth-dom_map(y)) * trunc(decaywidth-b+dom_map(y)); }
T u2_x(T x)      {  return  1. * trunc(a+decaywidth-dom_map(x)) * trunc(decaywidth-b+dom_map(x)); }
T u2_y(T y)      {  return  weight2*max(1.-exp(dom_map(y)),0.L) * trunc(a+decaywidth-dom_map(y))  * trunc(decaywidth-b+dom_map(y)); }

T f_t(T t)       {  return 0.; }
T f_x(T x)       {  return 0.; }
T f_y(T y)       {  return 0.; }


T
plot_initial_cond(const Basis2D& basis2d, const Coefficients<Lexicographical,T,Index2D> &u);

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
    T timestep_eps = 1e-6;
    //int maxiterations = 50;  T init_cgtol = 1e-2;
    int maxiterations =  1;  T init_cgtol = 1e-8;   // use maxiterations = 1 for "pure" sparse grid computation
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
    T a1 = 0.5*sigma1*sigma1/((b-a)*(b-a));
    T a2 = 0.5*sigma2*sigma2/((b-a)*(b-a));
    UniDirectionalLocalOpXOne2D  uniDirectionalOpXOne2D(localOp1D, a1);
    UniDirectionalLocalOpXTwo2D  uniDirectionalOpXTwo2D(localOp1D, a2);
    CompoundLocalOperator2D      localOp2D(uniDirectionalOpXOne2D,uniDirectionalOpXTwo2D);
    ThetaTimeStepLocalOperator2D localThetaTimeStepOp2D(theta,timestep,localOp2D);

    /// Initialization of preconditioner
    NoPreconditioner<T, Index2D> NoPrec;
    Preconditioner  Prec(basis2d, a1, a2, 1.);

    /// Initialization of integrals for initial condition and rhs
    DenseVectorT sing_pts_t, sing_pts_x(5), sing_pts_y(5);
    sing_pts_x = 0.1, 0.2, 0.3, 0.4, 0.5;
    sing_pts_y =  0.1, 0.2, 0.3, 0.4, 0.5;
    DenseMatrixT no_deltas, deltas_x, deltas_y;
    int order = 40;
    Function<T>                    fct_f_t(f_t,sing_pts_t);
    Function<T>                    fct_f_x(f_x,sing_pts_x), fct_f_y(f_y,sing_pts_y);
    Function<T>                    fct_u1_x(u1_x,sing_pts_x), fct_u2_x(u2_x,sing_pts_x);
    Function<T>                    fct_u1_y(u1_y,sing_pts_y), fct_u2_y(u2_y,sing_pts_y);
    RHSWithPeaks1D<T,PrimalBasis>  rhs_f_x(basis, fct_f_x, no_deltas, order);
    RHSWithPeaks1D<T,PrimalBasis>  rhs_f_y(basis, fct_f_y, no_deltas, order);
    RHSWithPeaks1D<T,PrimalBasis>  rhs_u1_x(basis, fct_u1_x, no_deltas, order);
    RHSWithPeaks1D<T,PrimalBasis>  rhs_u2_x(basis, fct_u2_x, no_deltas, order);
    RHSWithPeaks1D<T,PrimalBasis>  rhs_u1_y(basis, fct_u1_y, no_deltas, order);
    RHSWithPeaks1D<T,PrimalBasis>  rhs_u2_y(basis, fct_u2_y, no_deltas, order);

    Coefficients<Lexicographical,T,Index1D> rhs_f_x_data(SIZEHASHINDEX1D),
                                            rhs_f_y_data(SIZEHASHINDEX1D),
                                            rhs_u1_x_data(SIZEHASHINDEX1D), rhs_u1_y_data(SIZEHASHINDEX1D),
                                            rhs_u2_x_data(SIZEHASHINDEX1D), rhs_u2_y_data(SIZEHASHINDEX1D);

    AdaptiveSeparableRhsIntegral2D F_rhs(rhs_f_x, rhs_f_x_data, rhs_f_y, rhs_f_y_data);
    ThetaTimeStepRhs2d thetatimestep_F(fct_f_t,F_rhs);

    AdaptiveSeparableRhsIntegral2D F_u0_1(rhs_u1_x, rhs_u1_x_data, rhs_u1_y, rhs_u1_y_data);
    AdaptiveSeparableRhsIntegral2D F_u0_2(rhs_u2_x, rhs_u2_x_data, rhs_u2_y, rhs_u2_y_data);
    CompoundRhsIntegral2D          F_u0(F_u0_1, F_u0_2);



    /// Initialization of multi tree based adaptive wavelet Galerkin method
    ApproxL2AWGM2D approxL2_solver(basis2d, localThetaTimeStepOp2D, F_u0, NoPrec);
    approxL2_solver.setParameters(alpha, gamma, residualType, treeType, IsMW, false);

    ThetaTimeStepMultiTreeAWGM2D thetatimestep_solver(basis2d, localThetaTimeStepOp2D,
                                                          thetatimestep_F, Prec);
    thetatimestep_solver.setParameters(alpha, gamma, residualType, treeType, IsMW, false);

    ThetaSchemeMultiTreeAWGM2D thetascheme(thetatimestep_solver);

    ofstream convfile("conv_sol.dat");
    for (int j=7; j<=10; ++j) {
        size_t hms = 49157;
        Coefficients<Lexicographical,T,Index2D> u(hms);
        getSparseGridVector(basis2d, u, j, 0.L);

        //ofstream convfile("conv_u0.dat");
        //for (T eps = 1.; eps >= 1e-5; eps*=0.1) {
            //approxL2_solver.approxL2(u, eps, maxiterations);
            //plotScatterCoeff(u, basis2d, "coeff_u0", true);
            //T error = plot_initial_cond(basis2d, u);
            //convfile << u.size() << " " << error << endl;
        //}

        approxL2_solver.approxL2(u, 0.1*timestep_eps, maxiterations);
        //plot_initial_cond(basis2d, u);
        //plotScatterCoeff(u, basis2d, "coeff_u0", true);

        //for (int numOfTimesteps=1; numOfTimesteps<=64; numOfTimesteps*=2) {
        T max_error = 0.;
        for (int numOfTimesteps=256; numOfTimesteps<=256; numOfTimesteps*=2) {
            timestep = 1./numOfTimesteps;
            thetascheme.setParameters(theta, timestep, numOfTimesteps, timestep_eps, maxiterations,
                                      init_cgtol);
            thetascheme.solve(u);
            max_error = plot_current_solution(numOfTimesteps*timestep, basis2d, u);
        }
        convfile << u.size() << " " << j << " " << max_error << endl;
        cout << "sigma1  = " << sigma1 << ",  sigma2  = " << sigma2 << endl;
        cout << "weight1 = " << weight1 << ", weight2 = " << weight2 << endl;
        //plotScatterCoeff(u, basis2d, "coeff_sol", true);
    }



    return 0;
}

T
plot_initial_cond(const Basis2D& basis2d, const Coefficients<Lexicographical,T,Index2D> &u)
{
    stringstream plotfilename;
    plotfilename << "u0.dat";
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

            plotfile << x << " " << y << " " << u0(x,y) << " " << u1_x(x)*u1_y(y) + u2_x(x)*u2_y(y) << " " << val << endl;
            max_error = std::max(max_error, fabs(val-u0(x,y)));
        }
        plotfile << endl;
    }
    std::cerr << "Maximum error: " << max_error << std::endl;
    return max_error;
}

T
plot_current_solution(T time, const Basis2D& basis2d, const Coefficients<Lexicographical,T,Index2D> &u)
{
    Option1D<T,optiontype>  option1(optionparameters);
    Option1D<T,optiontype>  option2(optionparameters);


    //exact = std::exp(r*maturity)*option.value(processparameters, spot, 0);

    stringstream plotfilename;
    plotfilename << "theta_sol_" << time << ".dat";
    ofstream plotfile(plotfilename.str().c_str());

    T max_error = 0.;
    for (T x=0.; x<=1.; x+=0.01) {
        T spot1 = strike*std::exp(dom_map(x)+(0.5*sigma1*sigma1-r)*maturity);
        T exact1 = weight1 * std::exp(r*maturity)*option1.value(processparameters1, spot1, 0);
        for (T y=0.; y<=1.; y+=0.01) {
            T spot2 = strike*std::exp(dom_map(y)+(0.5*sigma2*sigma2-r)*maturity);
            T exact2 = weight2 * std::exp(r*maturity)*option2.value(processparameters2, spot2, 0);
            T val = 0.;
            for (const_coeff2d_it it=u.begin(); it!=u.end(); ++it) {
                XType xtype_x = (*it).first.index1.xtype, xtype_y = (*it).first.index2.xtype;
                int   j_x = (*it).first.index1.j, j_y = (*it).first.index2.j;
                long  k_x = (*it).first.index1.k, k_y = (*it).first.index2.k;
                val += (*it).second * basis2d.first.generator(xtype_x)(x,j_x,k_x,0) * basis2d.second.generator(xtype_y)(y,j_y,k_y,0);
                //std::cerr << "val = " << val << std::endl;
            }
            plotfile << x << " " << y << " " << exact1 + exact2 << " " << val << endl;
            if ((0.4 <= x && x <= 0.6) && (0.4 <= y && y <= 0.6)) {
                max_error = std::max(max_error, fabs(val-(exact1 + exact2)));
            }
        }
        plotfile << endl;
    }
    std::cerr << "Maximum error: " << max_error << std::endl;
    return max_error;
}

