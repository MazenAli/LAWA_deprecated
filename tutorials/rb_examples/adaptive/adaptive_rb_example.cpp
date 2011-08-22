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

using namespace std;
using namespace lawa;

typedef double T;


// Basis definitions
typedef Basis<T,Primal,Interval,Dijkema>                    IntervalBasis;
typedef TensorBasis2D<Adaptive,IntervalBasis,IntervalBasis> Basis2D;

// Operator Definitions
typedef WeightedHelmholtzOperator2D<T, Basis2D>                            WeightedLaplaceOp2D;
typedef H1NormPreconditioner2D<T, Basis2D>                                DiagPrec2D;
typedef WeightedAdaptiveHelmholtzOperator2D<T, Basis2D, DiagPrec2D>        WeightedAdaptHHOp2D;
typedef AdaptiveHelmholtzOperator2D<T, Basis2D, DiagPrec2D>                AdaptHHOp2D;
typedef RHS<T,Index2D, SeparableRHS2D<T, Basis2D>, DiagPrec2D>          AdaptRHS;

// Algorithm Definition
    // Class for calling the truth solver for snapshot calculations
typedef CompressionWeightedPDE2D<T, Basis2D>                              Compression;
typedef S_ADWAV_TruthSolver<T, Basis2D, Index2D, Compression>             S_AdwavSolver;
// Class containing all \calN-dependent data and functions
typedef AdaptiveRBTruth2D<T, Basis2D, S_AdwavSolver, Compression>         RBTruth;
    // Class containing only N-dependent data and functions
typedef RBModel2D<T, RBTruth>                                             RBModel;
    // The actual S_adwav solver, used as truth solver
typedef S_ADWAV<T, Index2D, Basis2D, RBTruth::Operator_LHS, RBTruth::Operator_RHS>    S_Adwav;

// Data type definitions
typedef Coefficients<Lexicographical,T,Index2D>                        Coeffs;
typedef flens::DenseVector<flens::Array<T> >                         DenseVectorT;
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
weight_Omega_x_1(T x)
{
    return (x <= 0.5) ? 1 : 0;
}

T
weight_Omega_x_2(T x)
{
    return (x >= 0.5) ? 1 : 0;
}

T
weight_Omega_y(T y)
{
    return 1;
}

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
weight_Forcing_x(T x)
{
    return (x >= 0.4) && (x <= 0.6) ? 1 : 0;    
}

T
weight_Forcing_y(T y)
{
    return (y >= 0.45) && (y <= 0.55) ? 1 : 0;    
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

    if (argc != 7) {
        cerr << "Usage " << argv[0] << " d d_ j0_x j0_y max_its mu1" << endl; 
        exit(1);
    }

    int d       = atoi(argv[1]);
    int d_      = atoi(argv[2]);
    int j0_x    = atoi(argv[3]);
    int j0_y    = atoi(argv[4]);
    int numIts  = atoi(argv[5]);
    T     mu1        = atof(argv[6]);
    
    /* Basis Initialization */
    IntervalBasis   basis_x(d, d_, j0_x);
    IntervalBasis   basis_y(d, d_, j0_y);
    basis_y.enforceBoundaryCondition<DirichletBC>();
    Basis2D         basis2d(basis_x, basis_y);
    
    DiagPrec2D             prec(basis2d);

    /* Model Initialization */
    RBModel rb_model;
    
        // We need a truth model, as we want to do truth solves
    RBTruth rb_truth(basis2d);
    rb_model.set_truthmodel(rb_truth);
    // Attach an inner product (here: scalar product in L2)
    AdaptHHOp2D h1norm(basis2d, 0., prec);
    rb_truth.attach_inner_product_op(h1norm);
    
        // Parameter vector
    std::vector<T> mu(1);
    mu[0] = mu1;
    rb_model.set_current_param(mu);
    std::vector<T> refmu(1);
    refmu[0] = 1.;
    rb_model.set_ref_param(refmu);
  
    /* Affine Decomposition Left Hand Side */
    DenseVectorT singpts_x(3), singpts_y(2);
    singpts_x = 0., 0.5, 1.;
    singpts_y = 0., 1.;
    Function<T>         w_x_1(weight_Omega_x_1, singpts_x);
    Function<T>         w_x_2(weight_Omega_x_2, singpts_x);
    Function<T>         w_y(weight_Omega_y, singpts_y);
        
    WeightedAdaptHHOp2D hh_1(basis2d, 0, w_x_1, w_y, prec, 1e-10);
    WeightedAdaptHHOp2D hh_2(basis2d, 0, w_x_2, w_y, prec, 1e-10);
    
    rb_model.truth->attach_A_q(theta_a_1, hh_1);
    rb_model.truth->attach_A_q(theta_a_2, hh_2);

    cout << "Q_a = " << rb_model.Q_a() << endl;
    cout << "Theta_a_1(mu) = " << (*rb_model.theta_a[0])(mu) << std::endl;

    /* Affine Decomposition Right Hand Side */
    DenseVectorT singpts_f_x(4), singpts_f_y(4);
    singpts_f_x = 0., 0.4, 0.6, 1.;
    singpts_f_y = 0., 0.45, 0.55, 1.;
    SeparableFunction2D<T> forcingFct(weight_Forcing_x, singpts_f_x, weight_Forcing_y, singpts_f_y);
    FullColMatrixT noDeltas;
    SeparableRHS2D<T, Basis2D> forcingIntegral(basis2d, forcingFct, noDeltas, noDeltas, 4);
    AdaptRHS F(forcingIntegral, prec);
    
    rb_model.truth->attach_F_q(theta_f_1, F);
    
    cout << "Q_f = " << rb_model.Q_f() << endl;
    
    /* Truth Solver: S_Adwav */
    
    T contraction   = 0.125;
    T threshTol     = 0.1;
    T cgTol         = 1e-6;
    T resTol        = 1e-4;
    
    SolverCall call = call_cg;
    S_AdwavSolver s_adwav_solver(rb_truth, call);
    s_adwav_solver.set_parameters(contraction, threshTol, cgTol, resTol, numIts);
    rb_truth.set_truthsolver(s_adwav_solver);
    
    /* Training (as yet: "manually") */
    
    // First Basis Function
    cout << "Solving for parameter mu = " << mu[0] << endl;
    Coeffs u = rb_model.truth->solver->truth_solve();
    rb_model.add_to_basis(u);
    
    cout << "N = " << rb_model.n_bf() << endl;
    cout << "F = " << rb_model.RB_F_vectors[0] << endl;
    cout << "A_1 = " << rb_model.RB_A_matrices[0] << endl;
    cout << "A_2 = " << rb_model.RB_A_matrices[1] << endl;
    cout << "    Snapshot 1: " << u.size() << " Basis Functions" << std::endl;
    ofstream snapshot1("Coeffs_Snapshot1.txt");
    snapshot1 << u << endl;
    snapshot1.close();
    plot2D(basis2d, u, prec, plot_dummy_fct, 0., 1., 0., 1., 0.02, "adaptive_snapshot_1");
    
    cout << "-----------------------" << endl << endl;
    
    // Second  Basis Function
    mu[0] = 0.5 * mu1;
    rb_model.set_current_param(mu);
    cout << "Solving for parameter mu = [" << mu[0] << ", " << mu[1] << "] "<< endl;
    u.clear();
    u = rb_model.truth->solver->truth_solve();
    rb_model.add_to_basis(u);
    
    cout << "N = " << rb_model.n_bf() << endl;
    cout << "F = " << rb_model.RB_F_vectors[0] << endl;
    cout << "A_1 = " << rb_model.RB_A_matrices[0] << endl;
    cout << "A_2 = " << rb_model.RB_A_matrices[1] << endl;
    cout << "    Snapshot 2: " << u.size() << " Basis Functions" << std::endl;
    ofstream snapshot2("Coeffs_Snapshot2.txt");
    snapshot2 << u << endl;
    snapshot2.close();
    plot2D(basis2d, u, prec, plot_dummy_fct, 0., 1., 0., 1., 0.02, "adaptive_snapshot_2");
    
    cout << "=======================" << endl << endl;

    // RB solve 
    mu[0] = 0.75 * mu1;
    rb_model.set_current_param(mu);
    cout << "Solving for parameter mu = " << mu[0] << endl;
    
    DenseVectorT rb_u = rb_model.RB_solve(rb_model.n_bf());
    cout << "    u_RB = " << rb_u << ", ResDualNorm = " << rb_model.residual_dual_norm(rb_u, mu)
                          <<  ", alpha_lb = " << rb_model.alpha_LB(mu) << endl;
    Coeffs u_N = rb_model.reconstruct_u_N(rb_u, rb_model.n_bf());
    
    cout << "    u_N          : " << u_N.size() << " Basis Functions" << std::endl;
    ofstream uN("Coeffs_u_N.txt");
    uN << u_N << endl;
    uN.close();
    plot2D(basis2d, u_N, prec, plot_dummy_fct, 0., 1., 0., 1., 0.02, "adaptive_u_N");


    return 0;
}
