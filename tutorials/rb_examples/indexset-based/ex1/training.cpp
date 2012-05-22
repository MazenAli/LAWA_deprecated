/*
 *  rb_example.cpp
 *  Xcode-LAWA
 *
 *  Created by Kristina Steih on 27.04.11.
 *  Copyright 2011 __MyCompanyName__. All rights reserved.
 *
 */

#include <iostream>
#include <iomanip>
#include <lawa/lawa.h>
#include "problem.h"

using namespace std;
using namespace lawa;

int main(int argc, char* argv[]) {

    /* PARAMETERS: */

    if (argc != 5) {
        cerr << "Usage " << argv[0] << " d d_ j0_x j0_y" << endl; 
        exit(1);
    }

    int d       = atoi(argv[1]);
    int d_      = atoi(argv[2]);
    int j0_x    = atoi(argv[3]);
    int j0_y    = atoi(argv[4]);
    
    /* Basis Initialization */
    IntervalBasis   basis_x(d, d_, j0_x);
    IntervalBasis   basis_y(d, d_, j0_y);
    basis_y.enforceBoundaryCondition<DirichletBC>();
    Basis2D         basis2d(basis_x, basis_y);
    
    DiagPrec2D             prec(basis2d);
    NoPrec2D               noprec;

    /* Model Initialization */
    RBModel rb_model;
    
    
        // We need a truth model, as we want to do truth solves
    bool use_inner_product_matrix = true;
    bool use_A_operator_matrices = false;
    RBTruth rb_truth(basis2d, prec, use_inner_product_matrix, use_A_operator_matrices);
    rb_model.set_truthmodel(rb_truth);
    
    // Attach an inner product (here: H1 semi norm)
    AdaptHHOp2D h1norm(basis2d, 1., noprec, 1e-10);
    rb_truth.attach_inner_product_op(h1norm);
    
        // Parameter vector
    std::vector<T> refmu(1);
    refmu[0] = reference_mu;
    rb_model.set_ref_param(refmu);
    rb_model.set_current_param(refmu);
  
    /* Affine Decomposition Left Hand Side */
    DenseVectorT singpts_x(3), singpts_y(2);
    singpts_x = 0., 0.5, 1.;
    singpts_y = 0., 1.;
    Function<T>         w_x_1(weight_Omega_x_1, singpts_x);
    Function<T>         w_x_2(weight_Omega_x_2, singpts_x);
    Function<T>         w_y(weight_Omega_y, singpts_y);
        
    WeightedAdaptHHOp2D hh_1(basis2d, 0, w_x_1, w_y, noprec);
    WeightedAdaptHHOp2D hh_2(basis2d, 0, w_x_2, w_y, noprec);
    
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
    AdaptRHS F(forcingIntegral, noprec);
    
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
    if(use_inner_product_matrix){
        rb_truth.assemble_inner_product_matrix(basisset);    
    }
    if(use_A_operator_matrices){
        rb_truth.assemble_A_operator_matrices(basisset);        
    }
    
    cout << "TEST - Solving for parameter mu = " << rb_model.get_current_param()[0] << endl;
    Coeffs u = rb_model.truth->truth_solve();
    plot2D(basis2d,u,noprec,plot_dummy_fct,0., 1., 0., 1., 0.01, "Test_Snapshot_ReferenceParameter");
    
    /* Training */
    
    std::vector<T> min_param, max_param, init_param;
    min_param.push_back(mu_min);
    max_param.push_back(mu_max);
    init_param.push_back(mu_init);
    std::vector<int> n_train_vec;
    n_train_vec.push_back(n_train);

    rb_model.set_min_param(min_param);
    rb_model.set_max_param(max_param);
    rb_model.generate_uniform_trainingset(n_train_vec);
    
    for(unsigned int i = 0; i < rb_model.Xi_train.size(); ++i) {
      for(unsigned int d = 0; d < rb_model.Xi_train[i].size(); ++d) {
          cout << rb_model.Xi_train[i][d] << endl;
      }
    }

    string trainingerrorfile = "Training.txt";
    
    rb_model.train_Greedy(min_param, Greedy_tol, Nmax, trainingerrorfile.c_str(), call);
    rb_model.write_RB_data();
    rb_model.write_basis_functions();
    rb_model.truth->write_riesz_representors();
    
    return 0;
}
