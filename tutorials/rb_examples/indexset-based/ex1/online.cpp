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
#include "problem.h"

using namespace std;
using namespace lawa;

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
