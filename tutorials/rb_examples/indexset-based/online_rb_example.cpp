/*
 *  rb_example.cpp
 *  Xcode-LAWA
 *
 *  Created by Kristina Steih on 27.04.11.
 *  Copyright 2011 __MyCompanyName__. All rights reserved.
 *
 */

#include <iostream>
#include <lawa/lawa.h>
#include <lawa/methods/rb/postprocessing/plotting.h>

using namespace std;
using namespace lawa;

typedef double T;


// Basis definitions
typedef Basis<T,Primal,Interval,Dijkema>                    IntervalBasis;
typedef TensorBasis2D<Adaptive,IntervalBasis,IntervalBasis> Basis2D;

// Operator Definitions
typedef WeightedHelmholtzOperator2D<T, Basis2D>                            WeightedLaplaceOp2D;
typedef H1NormPreconditioner2D<T, Basis2D>                                 DiagPrec2D;
typedef NoPreconditioner<T, Index2D>                                       NoPrec2D;
typedef WeightedAdaptiveHelmholtzOperator2D<T, Basis2D, NoPrec2D>          WeightedAdaptHHOp2D;
typedef AdaptiveHelmholtzOperator2D<T, Basis2D, NoPrec2D>                  AdaptHHOp2D;
typedef RHS<T,Index2D, SeparableRHS2D<T, Basis2D>, NoPrec2D>               AdaptRHS;

// Algorithm Definition
    // Class for calling the truth solver for snapshot calculations
typedef CompressionWeightedPDE2D<T, Basis2D>                              Compression;
typedef IndexsetTruthSolver<T, Basis2D, DiagPrec2D, 
                            Index2D, Compression>                         IndexsetSolver;
    // Class containing all \calN-dependent data and functions
typedef AdaptiveRBTruth2D<T, Basis2D, DiagPrec2D, 
                            IndexsetSolver, Compression>                  RBTruth;
    // Class containing only N-dependent data and functions
typedef RBModel2D<T, RBTruth>                                             RBModel;

// Data type definitions
typedef Coefficients<Lexicographical,T,Index2D>                       Coeffs;
typedef flens::DenseVector<flens::Array<T> >                          DenseVectorT;
typedef flens::GeMatrix<flens::FullStorage<T, cxxblas::ColMajor> >    FullColMatrixT;


/* Example: Thermal block  
 *
 *          - theta * u_xx = f on (0,1)^2
 *         u(x,0) = u(x,1) = 0
 *
 *        Here: theta = theta_1 on (0, 0.5) x (0,1)
 *              theta = theta_2 on (0,5, 1) x (0,1)
 *                  f = 1       on (0.4, 0.6) x (0.45, 0.55)
 */

T
theta_a_1(const std::vector<T>& mu)
{
    return mu[0];
}

T
theta_a_2(const std::vector<T>&)
{
    return 1;
}

T
theta_f_1(const std::vector<T>&)
{
    return 1.;
}

T
plot_dummy_fct(T x, T y)
{
    return 0;
}

int main(int argc, char* argv[]) {

    /* PARAMETERS: */

    if (argc != 8) {
        cerr << "Usage " << argv[0] << " d d_ j0_x j0_y mu N_test plot_h" << endl; 
        exit(1);
    }

    int d       = atoi(argv[1]);
    int d_      = atoi(argv[2]);
    int j0_x    = atoi(argv[3]);
    int j0_y    = atoi(argv[4]);
    T   mu      = atof(argv[5]);
    unsigned int N    = atoi(argv[6]);
    T h         = atof(argv[7]);
    
    /* Basis Initialization */
    IntervalBasis   basis_x(d, d_, j0_x);
    IntervalBasis   basis_y(d, d_, j0_y);
    basis_y.enforceBoundaryCondition<DirichletBC>();
    Basis2D         basis2d(basis_x, basis_y);
    
    DiagPrec2D             prec(basis2d);
    NoPrec2D               noprec;

    /* Model Initialization */
    RBModel rb_model;
    
        // Parameter vector
    std::vector<T> refmu(1);
    refmu[0] = 1.;
    std::vector<T> curr_mu(1);
    curr_mu[0] = mu;
    rb_model.set_ref_param(refmu);
    rb_model.set_current_param(curr_mu);
  
    /* Attach Theta-Functions */
    
    rb_model.attach_theta_a_q(theta_a_1);
    rb_model.attach_theta_a_q(theta_a_2);

    cout << "Q_a = " << rb_model.Q_a() << endl;
    
    rb_model.attach_theta_f_q(theta_f_1);
    
    cout << "Q_f = " << rb_model.Q_f() << endl;
    
    /* Read RB_data that has been computed and stored in offline phase */
    
    rb_model.read_RB_data();
    rb_model.read_basis_functions();

    /* RB_solve and plot */
    
    DenseVectorT u_RB = rb_model.RB_solve(N, curr_mu);
    cout << "ErrorEst at " << curr_mu[0] << ": " << rb_model.residual_dual_norm(u_RB, curr_mu) / rb_model.alpha_LB(curr_mu) << endl;
    
    Coeffs u_approx = rb_model.reconstruct_u_N(u_RB, N);
    
    stringstream filename, filename_plot;
    filename_plot << "u_N_" << N << "_mu_" << mu << "_plot";
    filename << "u_N_" << N << "_mu_" << mu << ".dat";
    cout << "Save coefficients ... " << endl;
    saveCoeffVector2D(u_approx, basis2d, filename.str().c_str());
    cout << "    ... done " << endl;
    cout << "Plotting function ... " << endl;
    plot2D<T>(basis2d, u_approx, noprec, plot_dummy_fct, 0, 1, 0, 1, h, filename_plot.str().c_str());
    cout << "    ... done " << endl;
    
    return 0;
}
