#include <fstream>
#include <iostream>
#include <flens/flens.h>
#include <lawa/lawa.h>
#include <applications/finance/fourierpricing/fourierpricing.h>
#include <applications/finance/initialconditions/initialcondition1d.h>
#include <applications/finance/options/options.h>
#include <applications/finance/operators/operators.h>
#include <applications/finance/processes/processes.h>
#include <applications/finance/righthandsides/righthandsides.h>

using namespace std;
using namespace lawa;

typedef double T;
typedef flens::GeMatrix<flens::FullStorage<T, cxxblas::ColMajor> >  FullColMatrixT;
typedef flens::SparseGeMatrix<flens::CRS<T,flens::CRS_General> >    SparseMatrixT;
typedef flens::DiagonalMatrix<T>                                    DiagonalMatrixT;
typedef flens::DenseVector<flens::Array<T> >                        DenseVectorT;

const OptionType1D optiontype = Put;
OptionParameters1D<T,Put> optionparameters(80.0, 1., false);

const ProcessType1D  processtype  = BlackScholes;
ProcessParameters1D<T,BlackScholes>   processparameters(0.04, 0.2);
//ProcessParameters1D<T,CGMY>             processparameters(0.04, 1., 2.4, 4.5, 1.8 );
//ProcessParameters1D<T,CGMYe>          processparameters(0.04, 1., 2.4, 4.5, 1.8, 0.1 );


typedef Basis<T,Primal,Interval,Dijkema>        Basis1D;

typedef WeightedL2ScalarProduct1D<T, Basis1D>               ScalarProduct1D;
typedef DiagonalMatrixPreconditioner1D<T, Basis1D,
                                       ScalarProduct1D>     WeightedL2Preconditioner1D;
typedef FinanceOperator1D<T, processtype, Basis1D>          FinanceOp;

typedef OptionRHS1D<T, optiontype, processtype, Basis1D>    OptionRhs;

typedef PayoffInitialCondition1D<optiontype, Basis1D>       PayoffInitCond;

// TimeStepping Methods
typedef ThetaScheme1D_LTI<T, Basis1D, FinanceOp, OptionRhs> Theta_Use_Excess_To_Payoff;
typedef TimeStepping<T, Theta_Use_Excess_To_Payoff>         TimeStepper_Use_Excess_To_Payoff;

typedef ThetaScheme1D_LTI<T, Basis1D, FinanceOp,
                          HomogeneousRHS<T> >               Theta_Use_Init_Cond;
typedef TimeStepping<T, Theta_Use_Init_Cond>                TimeStepper_Use_Init_Cond;


template<typename T, OptionType1D OType, ProcessType1D PType, typename Basis>
void
ComputeL2ErrorAndPlotSolution(Option1D<T,OType> &option,
                              ProcessParameters1D<T,PType> &processparameters,
                              const Basis &basis, int J,
                              const DenseVectorT &u0, const DenseVectorT &u,
                              T &L2error, T &Linftyerror, T eta, T R1, T R2,
                              int use_excess_to_payoff);

int
main(int argc, char *argv[])
{
    cout.precision(8);
    if (argc != 6) {
        cerr << "usage: " << argv[0] << " j0 J theta timesteps use_excess_to_payoff" << endl;
        exit(1);
    }
    int d=2, d_=2;
    int j0                   = atoi(argv[1]);
    int j_max                = atoi(argv[2]);
    T theta                  = T(atof(argv[3]));
    int timesteps            = atoi(argv[4]);
    int use_excess_to_payoff = atoi(argv[5]);

    if (theta < 0.5) {
        cout << "theta should be larger than 0.5!" << endl;
        exit(1);
    }

    T                       eta=2.;
    int                     order=20;

    T                       timestep = optionparameters.maturity/T(timesteps);
    T                       R1=4.;
    T                       R2=4.;
    Basis1D                 basis(d,d_,j0);
    basis.enforceBoundaryCondition<DirichletBC>();

    WeightedL2ScalarProduct1D<T, Basis1D> wL2scalarproduct(basis, 2., 0., 1., 8);

    Option1D<T,optiontype>  option(optionparameters);

    for (int J=j0; J<=j_max; ++J) {
        Timer time;
        time.start();
        DenseVectorT u(basis.mra.rangeI(J)), u0(basis.mra.rangeI(J));
        if (use_excess_to_payoff==1) {
            FinanceOp                         finance_op(basis, processparameters, eta, R1, R2, order, J);
            OptionRhs                         rhs(optionparameters, processparameters, basis,
                                                  R1, R2);
            Theta_Use_Excess_To_Payoff        scheme(theta, basis, finance_op, rhs, true,
                                                     false, 0., 1e-12, eta, R1, R2, order);
            TimeStepper_Use_Excess_To_Payoff  timestepmethod(scheme, timestep, timesteps, J);
            u = timestepmethod.solve(u0, false);
        }
        else {
            ScalarProduct1D            scalarproduct(basis, eta, R1, R2, order);
            WeightedL2Preconditioner1D preconditioner(scalarproduct);
            PayoffInitCond             payoff_initcond(option,basis,eta,R1,R2);
            Assembler1D<T, Basis1D>    assembler(basis);
            SparseMatrixT   M      =   assembler.assembleStiffnessMatrix(scalarproduct, J);
            DiagonalMatrixT P      =   assembler.assemblePreconditioner(preconditioner, J);
            DenseVectorT    rhs_u0 =   assembler.assembleRHS(payoff_initcond, J);

            cout << pcg(P, M, u0, rhs_u0) << " pcg iterations required." << endl;

            FinanceOp                 finance_op(basis, processparameters, eta, R1, R2, order, J);
            HomogeneousRHS<T>         homogeneousrhs;
            Theta_Use_Init_Cond       scheme(theta, basis, finance_op, homogeneousrhs, true, false,
                                             0., 1e-12, eta, R1, R2, order);

            TimeStepper_Use_Init_Cond timestepmethod(scheme, timestep, timesteps, J);
            u = timestepmethod.solve(u0, false);

        }
        time.stop();
        T comp_time = time.elapsed();

        T L2error = 0.;
        T Linftyerror = 0.;
        ComputeL2ErrorAndPlotSolution(option, processparameters, basis, J, u0, u,
                                      L2error, Linftyerror, eta, R1, R2, use_excess_to_payoff);
        cout      << J << " " << basis.mra.cardI(J) << " " << comp_time << " "
                  << L2error << " " << Linftyerror << endl;
    }



    return 0;
}

template<typename T, OptionType1D OType, ProcessType1D PType, typename Basis>
void
ComputeL2ErrorAndPlotSolution(Option1D<T,OType> &option,
                              ProcessParameters1D<T,PType> &processparameters,
                              const Basis &basis, int J,
                              const DenseVectorT &u0, const DenseVectorT &u,
                              T &L2error, T &Linftyerror, T eta, T R1, T R2, int use_excess_to_payoff)
{
    std::stringstream filename;
    filename << "tmp.txt";
    std::ofstream plotFile(filename.str().c_str());
    std::stringstream filename2;
    filename2 << "tmp2.txt";
    std::ofstream plotFile2(filename2.str().c_str());
    T maturity = option.optionparameters.maturity;
    T strike   = option.optionparameters.strike;
    T r        = processparameters.r;

    L2error     = 0.;
    Linftyerror = 0.;
    T h = pow2i<T>(-J-2)*(R1+R2);
    T delta = pow2i<T>(0);
    for (T x=-delta*R1; x<=delta*R2; x+=h) {
        T P_u0   = (1./sqrt(R1+R2))*evaluate(basis, J, u0, (x+R1)/(R1+R2), 0);
        T approx = (1./sqrt(R1+R2))*evaluate(basis, J, u, (x+R1)/(R1+R2), 0);
        T exact = 0.;
        T spot = strike*std::exp(x-r*maturity);
        exact = std::exp(r*maturity)*option.value(processparameters, spot, 0);

        if (use_excess_to_payoff==1) exact -= option.payoff(strike*exp(x));

        if ((fabs(x+delta*R1)<1e-12) || (fabs(x-delta*R2) < 1e-12)) {
            L2error += 0.5*h*std::pow(approx-exact,2.);
        }
        else    {
            L2error += h*std::pow(approx-exact,2.);
        }
        Linftyerror = std::max(Linftyerror, fabs(exact-approx));

        plotFile  << x << " " << approx << " " << exact << " "
                  << P_u0 << " " << option.payoff_log(x) << endl;

        T weight = exp(-eta*fabs(x));
        plotFile2 << x << " " << weight*approx << " " << weight*exact << " "
                  << weight*P_u0 << " " << weight*option.payoff_log(x) << endl;
    }
    plotFile.close();
    plotFile2.close();
    L2error = sqrt(L2error);
}



