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
typedef H1NormPreconditioner2D<T, Basis2D>                                DiagPrec2D;
typedef WeightedAdaptiveHelmholtzOperator2D<T, Basis2D, DiagPrec2D>        WeightedAdaptHHOp2D;
typedef AdaptiveHelmholtzOperator2D<T, Basis2D, DiagPrec2D>                AdaptHHOp2D;
typedef RHS<T,Index2D, SeparableRHS2D<T, Basis2D>, DiagPrec2D>             AdaptRHS;

// Algorithm Definition
    // Class for calling the truth solver for snapshot calculations
typedef CompressionWeightedPDE2D<T, Basis2D>                              Compression;
typedef IndexsetTruthSolver<T, Basis2D, Index2D, Compression>             IndexsetSolver;
    // Class containing all \calN-dependent data and functions
typedef AdaptiveRBTruth2D<T, Basis2D, IndexsetSolver, Compression>        RBTruth;
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

    if (argc != 6) {
        cerr << "Usage " << argv[0] << " d d_ j0_x j0_y max_its" << endl; 
        exit(1);
    }

    int d       = atoi(argv[1]);
    int d_      = atoi(argv[2]);
    int j0_x    = atoi(argv[3]);
    int j0_y    = atoi(argv[4]);
    int numIts  = atoi(argv[5]);
    
    /* Basis Initialization */
    IntervalBasis   basis_x(d, d_, j0_x);
    IntervalBasis   basis_y(d, d_, j0_y);
    basis_y.enforceBoundaryCondition<DirichletBC>();
    Basis2D         basis2d(basis_x, basis_y);
    
    DiagPrec2D             prec(basis2d);

    /* Model Initialization */
    RBModel rb_model;
    
        // Attach an inner product (here: H1 semi norm)
    AdaptHHOp2D h1norm(basis2d, 1., prec, 1e-10);
    rb_model.attach_inner_product_op(h1norm);
    
        // We need a truth model, as we want to do truth solves
    RBTruth rb_truth(basis2d);
    rb_model.set_truthmodel(rb_truth);
    
        // Parameter vector
    std::vector<T> refmu(1);
    refmu[0] = 1.;
    rb_model.set_ref_param(refmu);
    rb_model.set_current_param(refmu);
  
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
    
    /* Truth Solver: on fixed indexset */
    IndexSet<Index2D> basisset;
    
    stringstream filename;
    filename << "IndexSet.dat";
    readIndexSet2D<T>(basisset, filename.str().c_str(), true);      
    std::cout << "Read file " << filename.str() << std::endl;
    std::cout << "Basis Set: " << basisset.size() << " functions" << std::endl;

    SolverCall call = call_cg;
    IndexsetSolver indexset_solver(basisset, rb_truth, call);
    rb_truth.set_truthsolver(indexset_solver);
    
    /* Training */
    std::vector<T> min_param, max_param;
    min_param.push_back(0.1);
    max_param.push_back(2);
    std::vector<int> n_train;
    n_train.push_back(30);
    rb_model.set_min_param(min_param);
    rb_model.set_max_param(max_param);
    rb_model.generate_uniform_trainingset(n_train);
    
    for(unsigned int i = 0; i < rb_model.Xi_train.size(); ++i) {
      for(unsigned int d = 0; d < rb_model.Xi_train[i].size(); ++d) {
          cout << rb_model.Xi_train[i][d] << endl;
      }
    }
    
    T tol = 1e-10;
    int Nmax = 7;
    rb_model.train_Greedy(min_param, tol, Nmax);
    rb_model.write_RB_data();
    rb_model.write_basis_functions();
    rb_model.truth->write_riesz_representors();
    
    return 0;
}
