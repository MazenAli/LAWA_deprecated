#include <iostream>
#include <lawa/lawa.h>

#include <applications/finance/initialconditions/initialconditions.h>
#include <applications/finance/options/options.h>
#include <applications/finance/operators/operators.h>
#include <applications/finance/processes/processes.h>

using namespace std;
using namespace lawa;

typedef long double T;

typedef flens::GeMatrix<flens::FullStorage<T, cxxblas::ColMajor> >  DenseMatrixT;
typedef flens::DenseVector<flens::Array<T> >                        DenseVectorT;

typedef Basis<T,Orthogonal,Interval,Multi>                          PrimalBasis;
typedef TensorBasis2D<Adaptive,PrimalBasis,PrimalBasis>             Basis2D;

/* ************************************ */
/* *** Typedefs for financial model *** */
/* ************************************ */


T strike = 1.;
T maturity = 1.;
T weight1 = 0.5, weight2 = 0.5;

/*
const OptionTypenD optiontype = BasketPut;
OptionParameters2D<T,BasketPut> optionparameters(strike, maturity, weight1, weight2, false);
typedef PayoffIntegral2D<FullGridGL,Basis2D,TruncatedBasketPutOption2D<T> > PayoffIntegral;
*/

const OptionTypenD optiontype = SumOfPuts;
OptionParameters2D<T,SumOfPuts> optionparameters(strike, strike, maturity, weight1, weight2, false);
typedef PayoffIntegral2D<FullGridGL,Basis2D,TruncatedSumOfPutsOption2D<T> > PayoffIntegral;

const ProcessType2D  processtype  = CGMYeUnivariateJump2D;
//T r = 0.04; T sigma1 = 0.3, sigma2 = 0.2, rho = 0.;
//T u11 = 1., u12 = 0., u21 = 0., u22 = 1.;
T r = 0.;
T sigma1 = 0.3;
T sigma2 = 0.2;
T rho = 0.;
T k_C1 = 1., k_G1 = 5., k_M1 = 5., k_Y1 = 0.8;
T k_C2 = 1., k_G2 = 5., k_M2 = 5., k_Y2 = 0.8;


T u11 = 1., u12 = 0., u21 = 0., u22 = 1.;
T s1  = sigma1*sigma1, s2  = sigma2*sigma2;
T    critical_line_x1 = 0.6;
bool critical_above_x1 = true;


ProcessParameters2D<T,CGMYeUnivariateJump2D>   processparameters(r, sigma1, sigma2, rho,
                                                                 k_C1,  k_G1, k_M1, k_Y1,
                                                                 k_C2,  k_G2, k_M2, k_Y2);

/* ********************************************* */
/* *** Typedefs for numerical discretization *** */
/* ********************************************* */

//typedef OptimizedH1Preconditioner2D<T,Basis2D>                      Preconditioner;
typedef FinanceOperator2D<CGMYeUnivariateJump2D, Basis2D>           CGMYeOp2D;
typedef ThetaTimeStepLocalOperator<Index2D,CGMYeOp2D>               ThetaTimeStepLocalOperator2D;
typedef DiagonalMatrixPreconditioner2D<T,Basis2D,CGMYeOp2D>         Preconditioner;

//Righthandsides definitions (separable)
typedef RHSWithPeaks1D<T,PrimalBasis>                               Rhs1D;
typedef AdaptiveSeparableRhs<T,Index2D,Rhs1D,Rhs1D >                AdaptiveSeparableRhsIntegral2D;
typedef ThetaTimeStepSeparableRHS<T,Index2D,
                                  AdaptiveSeparableRhsIntegral2D,
                                  ThetaTimeStepLocalOperator2D>     ThetaTimeStepRhs2d;
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


T f_t(T t)       {  return 0.; }
T f_x(T x)       {  return 0.; }
T f_y(T y)       {  return 0.; }

T
computeLinftyError(const Basis2D &basis2d, T left_x1, T right_x1, T left_x2, T right_x2,
                   const Coefficients<Lexicographical,T,Index2D> &u,T delta, int j,
                   Option2D<T,optiontype> &option2d,
                   ProcessParameters2D<T,CGMYeUnivariateJump2D> &processparameters);

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
    T gamma = 0.025;
    const char* residualType = "standard";
    const char* treeType = "sparsetree"; //"gradedtree";
    bool IsMW = true;
    size_t hashMapSize = 196613;
    T R1_1 = 4., R2_1 = 4.;
    T left_x1 = -R1_1, right_x1 = R2_1;
    T R1_2 = 4., R2_2 = 4.;
    T left_x2 = -R1_2, right_x2 = R2_2;
    T delta = 0.05;

    T theta = 0.5;
    T timestep_eps = 1e-6;
    int maxiterations =  1;  T init_cgtol = 1e-9;   // use maxiterations = 1 for "pure" sparse grid computation
    int numOfTimesteps = 16;
    T timestep = maturity/numOfTimesteps;

    int numOfMCRuns = 100000;

    int order = 4;

    bool useRefPrices = false;

    Timer time;

    /// Basis initialization
    PrimalBasis     basis(d,j0);
    basis.enforceBoundaryCondition<DirichletBC>();
    Basis2D         basis2d(basis,basis);

    cout << "Process parameters: " << processparameters << endl;
    CGMYeOp2D                    cgmyeOp2D(basis2d, processparameters,
                                           R1_1, R2_1, R1_2, R2_2, 10);
    ThetaTimeStepLocalOperator2D localThetaTimeStepOp2D(theta,timestep,cgmyeOp2D);

    /// Initialization of preconditioner
    NoPreconditioner<T, Index2D> NoPrec;
    //Preconditioner  Prec(basis2d, sigma1*sigma1, sigma2*sigma2, 1.);
    Preconditioner  Prec(cgmyeOp2D);


    /// Initialization of integrals for initial condition and rhs
    DenseVectorT sing_pts_t, sing_pts_x(5), sing_pts_y(5);
    sing_pts_x = 0.1, 0.2, 0.3, 0.4, 0.5;
    sing_pts_y =  0.1, 0.2, 0.3, 0.4, 0.5;
    DenseMatrixT no_deltas, deltas_x, deltas_y;
    Function<T>                    fct_f_t(f_t,sing_pts_t);
    Function<T>                    fct_f_x(f_x,sing_pts_x), fct_f_y(f_y,sing_pts_y);
    RHSWithPeaks1D<T,PrimalBasis>  rhs_f_x(basis, fct_f_x, no_deltas, order);
    RHSWithPeaks1D<T,PrimalBasis>  rhs_f_y(basis, fct_f_y, no_deltas, order);

    Coefficients<Lexicographical,T,Index1D> rhs_f_x_data(SIZEHASHINDEX1D),
                                            rhs_f_y_data(SIZEHASHINDEX1D);

    AdaptiveSeparableRhsIntegral2D F_rhs(rhs_f_x, rhs_f_x_data, rhs_f_y, rhs_f_y_data);
    ThetaTimeStepRhs2d thetatimestep_F(fct_f_t,F_rhs,localThetaTimeStepOp2D);

    /// Initialization of integrals for initial condition and rhs
    Option2D<T,optiontype>         option2d(optionparameters);

    //TruncatedBasketPutOption2D<T> truncatedoption2d;
    TruncatedSumOfPutsOption2D<T> truncatedoption2d;
    truncatedoption2d.setOption(option2d);
    truncatedoption2d.setTransformation(u11, u21, u12, u22);
    truncatedoption2d.setTruncation(left_x1, right_x1, left_x2, right_x2, 0, 0.1, 100.);
    truncatedoption2d.setCriticalLine_x1(critical_line_x1, critical_above_x1);

    PayoffIntegral payoffIntegral(basis2d, truncatedoption2d,
                                  left_x1, right_x1, left_x2, right_x2, true, 0.05, order);

    Coefficients<Lexicographical,T,Index2D> u(hashMapSize), f(hashMapSize);


    for (int j=0; j<=J; ++j) {

        getSparseGridVector(basis2d, u, j, (T)0.);

        cerr << "Computation of initial condition started." << endl;
        time.start();
        int count = 0;
        for (coeff2d_it it=u.begin(); it!=u.end(); ++it) {
            coeff2d_it it_f = f.find((*it).first);
            if (it_f != f.end()) {
                (*it).second = (*it_f).second;
            }
            else {
                T tmp = payoffIntegral((*it).first);
                f[(*it).first] = tmp;
                (*it).second = tmp;
            }
            ++count;
            if (count%100==0) cout << "count: " << count << " / " << u.size() << endl;

        }
        timestep = maturity/numOfTimesteps;

        /// Initialization of multi tree based adaptive wavelet Galerkin method
        ThetaTimeStepMultiTreeAWGM2D thetatimestep_solver(basis2d, localThetaTimeStepOp2D,
                                                              thetatimestep_F, Prec);
        thetatimestep_solver.setParameters(alpha, gamma, residualType, treeType, IsMW, false,
                                           hashMapSize);

        ThetaSchemeMultiTreeAWGM2D thetascheme(thetatimestep_solver);
        thetascheme.setParameters(theta, timestep, numOfTimesteps, timestep_eps, maxiterations,
                                  init_cgtol, 0);
        int avDof = 0, maxDof = 0., terminalDof;
        thetascheme.solve(u, avDof, maxDof, terminalDof, j);
        cerr << "Computation of u has finished." << endl;
        T maxerror = 0., maxerror1 = 0., maxerror2 = 0.;
        computeLinftyError(basis2d, left_x1, right_x1, left_x2, right_x2, u, 1., j, option2d,
                           processparameters);
        cerr << "Computation of errors has finished." << endl;
        //computeReferencePrice(basis2d, left_x1, right_x1, left_x2, right_x2,
        //                      -0.1, 0.1, -0.1, 0.1, 0.02, 0.02, u, j, option2d, processparameters);
        cerr << "Computation of reference prices has finished." << endl;

    }

    return 0;
}

T
evaluate(const Basis2D &basis2d, T left_x1, T right_x1, T left_x2, T right_x2,
         const Coefficients<Lexicographical,T,Index2D> &v, T x1, T x2)
{
    T RightmLeft_x1 = right_x1-left_x1, SqrtRightmLeft_x1 = std::sqrt(right_x1-left_x1);
    T RightmLeft_x2 = right_x2-left_x2, SqrtRightmLeft_x2 = std::sqrt(right_x2-left_x2);

    T ret = 0.;
    for (const_coeff2d_it it=v.begin(); it!=v.end(); ++it) {
        int   j1 = (*it).first.index1.j,     j2 = (*it).first.index2.j;
        int   k1 = (*it).first.index1.k,     k2 = (*it).first.index2.k;
        XType e1 = (*it).first.index1.xtype, e2 = (*it).first.index2.xtype;

        T val_x1 = (1./SqrtRightmLeft_x1) * basis2d.first.generator(e1).operator()((x1-left_x1)/(RightmLeft_x1),j1,k1,0);
        T val_x2 = (1./SqrtRightmLeft_x2) * basis2d.second.generator(e2).operator()((x2-left_x2)/(RightmLeft_x2),j2,k2,0);

        ret += (*it).second * val_x1 * val_x2;
    }
    return ret;
}

T
computeLinftyError(const Basis2D &basis2d, T left_x1, T right_x1, T left_x2, T right_x2,
                   const Coefficients<Lexicographical,T,Index2D> &u,T delta, int j,
                   Option2D<T,optiontype> &option2d,
                   ProcessParameters2D<T,CGMYeUnivariateJump2D> &processparameters)
{
    std::stringstream filename;
    if (optiontype == BasketPut) {
        filename << "cgmye2d_basketput_sg_" << j
                 << "_strike_" << strike << "_maturity_" << maturity << "_weight1_" << weight1 << "_weight2_" << weight2
                 << processparameters << ".txt";
    }
    else if (optiontype == SumOfPuts) {
        filename << "cgmye2d_sumofputs_sg_" << j
                 << "_strike_" << strike << "_maturity_" << maturity << "_weight1_" << weight1 << "_weight2_" << weight2
                 << processparameters << ".txt";
    }
    else {
        std::cerr << "Unknown option type" << std::endl; exit(1);
    }
    std::ofstream plotfile(filename.str().c_str());
    plotfile.precision(16);

    T maxerror = 0.;
    T h1 = (delta*right_x1-delta*left_x1)/25.;
    T h2 = (delta*right_x2-delta*left_x2)/25.;
    //for (T x1=left_x1; x1<=right_x1; x1+=0.03125) {
    for (T x1=delta*left_x1; x1<=delta*right_x1; x1+=h1) {
        //for (T x2=left_x2; x2<=right_x2; x2+=0.03125) {
        for (T x2=delta*left_x2; x2<=delta*right_x2; x2+=h2) {
            T S1    = strike*std::exp(x1);
            T S2    = strike*std::exp(x2);
            T exact = option2d.value(processparameters,S1,S2,0);
            T x1hat = u11*x1+u12*x2-u11*(0.5*sigma1*sigma1-r)*maturity-u12*(0.5*sigma2*sigma2-r)*maturity;
            T x2hat = u21*x1+u22*x2-u21*(0.5*sigma1*sigma1-r)*maturity-u22*(0.5*sigma2*sigma2-r)*maturity;
            T payoff = option2d.payoff(strike*exp(x1),strike*exp(x2));
            T approx =std::exp(-r*maturity)*evaluate(basis2d, left_x1, right_x1, left_x2, right_x2, u, x1hat, x2hat);
            plotfile << x1 << " " << x2 << " " << exact << " " << approx << " " << payoff << endl;
            maxerror = std::max(maxerror, fabs(approx-exact));
        }
        plotfile << endl;
    }
    plotfile.close();

    return maxerror;
}


