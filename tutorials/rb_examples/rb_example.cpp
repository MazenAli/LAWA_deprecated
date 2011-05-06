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
typedef TensorBasis2D<Uniform,IntervalBasis,IntervalBasis>  Basis2D;

// Operator definitions
typedef H1NormPreconditioner2D<T, Basis2D>							DiagPrec2D;
typedef UniformTruthSolver2D<T, Basis2D, DiagPrec2D>				UniformTruthSolver;
typedef UniformRBTruth2D<T, UniformTruthSolver>                     RBTruth;
typedef RBModel2D<T, RBTruth>                                       RBModel;

// FLENS definitions
typedef flens::DenseVector<flens::Array<T> > 						DenseVectorT;
typedef flens::GeMatrix<flens::FullStorage<T, cxxblas::ColMajor> >	FullColMatrixT;

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
theta_a_1(std::vector<T>& mu)
{
	return mu[0];
}

T
theta_a_2(std::vector<T>&)
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
theta_f_1(std::vector<T>&)
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
		cerr << "Usage " << argv[0] << " d d_ j0_x j0_y J_x J_y mu1" << endl; 
        exit(1);
	}

	int d       = atoi(argv[1]);
    int d_      = atoi(argv[2]);
	int j0_x    = atoi(argv[3]);
	int j0_y    = atoi(argv[4]);
	int J_x	    = atoi(argv[5]);
	int J_y	    = atoi(argv[6]);
	T mu1	    = atof(argv[7]);
    
    /* Basis Initialization */
    
	IntervalBasis   basis_x(d, d_, j0_x);
    IntervalBasis   basis_y(d, d_, j0_y);
	basis_y.enforceBoundaryCondition<DirichletBC>();
	Basis2D         basis2d(basis_x, basis_y);
    
	DiagPrec2D 			prec(basis2d);
    
    RBModel rb_model;
    RBTruth rb_truth;
    UniformTruthSolver truth_solver(basis2d, J_x, J_y, prec);
    rb_truth.set_truthsolver(truth_solver);
    rb_model.set_truthmodel(rb_truth);
    
    std::vector<T> mu(1);
    mu[0] = mu1;
    rb_model.set_current_param(mu);
    
    /* Affine Decomposition Left Hand Side */
    
    DenseVectorT singpts_x(3), singpts_y(2);
    singpts_x = 0., 0.5, 1.;
    singpts_y = 0., 1.;
    Function<T> w_x_1(weight_Omega_x_1, singpts_x);
    Function<T> w_x_2(weight_Omega_x_2, singpts_x);
    Function<T> w_y(weight_Omega_y, singpts_y);
    WeightedLaplaceOperator2D<T, Basis2D> hh_1(basis2d, w_x_1, w_y);
    WeightedLaplaceOperator2D<T, Basis2D> hh_2(basis2d, w_x_2, w_y);
    
    rb_truth.attach_A_q(theta_a_1, hh_1);
    rb_truth.attach_A_q(theta_a_2, hh_2);

	cout << "Q_a = " << rb_model.Q_a() << endl;
    cout << "Theta_a_1(mu) = " << (*rb_model.theta_a[0])(mu) << std::endl;
    
    /* Affine Decomposition Right Hand Side */
    
    DenseVectorT singpts_f_x(4), singpts_f_y(4);
    singpts_f_x = 0., 0.4, 0.6, 1.;
    singpts_f_y = 0., 0.45, 0.55, 1.;
    SeparableFunction2D<T> forcingFct(weight_Forcing_x, singpts_f_x, weight_Forcing_y, singpts_f_y);
    FullColMatrixT noDeltas;
    SeparableRHS2D<T, Basis2D> F(basis2d, forcingFct, noDeltas, noDeltas, 4);
    
    rb_truth.attach_F_q(theta_f_1, F);
    
    cout << "Q_f = " << rb_model.Q_f() << endl;
    
    /* Solution */
    
    Timer timer;
    timer.start();
    Coefficients<Lexicographical,T,Index2D> u = rb_model.truth->solver->truth_solve();
    timer.stop();
    cout << "Time: " << timer.elapsed() << " seconds " << endl;
    cout << u << endl;
    
    NoPreconditioner<T, Index2D> noprec;
    plot2D(basis2d, u, noprec, plot_dummy_fct, 0., 1., 0., 1., 0.01, "rb_example.txt");
    
	return 0;
}