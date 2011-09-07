/*
 * indexset_rb_output_test.tcc
 *
 *  Created on: 28.07.2011
 *      Author: s_sglas
 */
//Erweiterung (und Veränderung) des indexset_rb_test.tcc

#include <iostream>
#include <lawa/lawa.h>
#include <lawa/methods/rb/postprocessing/plotting.h>
#include <lawa/methods/rb/outputs/averageoutput.h>

using namespace std;
using namespace lawa;

typedef double T;


// Basis definitions
typedef Basis<T,Primal,Interval,Dijkema>                    IntervalBasis;
typedef TensorBasis2D<Adaptive,IntervalBasis,IntervalBasis> Basis2D;

// Operator Definitions
typedef WeightedHelmholtzOperator2D<T, Basis2D>                          WeightedLaplaceOp2D;
typedef H1NormPreconditioner2D<T, Basis2D>                               DiagPrec2D;
typedef NoPreconditioner<T, Index2D>                                     NoPrec2D;
typedef WeightedAdaptiveHelmholtzOperator2D<T, Basis2D, NoPrec2D>        WeightedAdaptHHOp2D;
typedef AdaptiveHelmholtzOperator2D<T, Basis2D, NoPrec2D>                AdaptHHOp2D;
typedef RHS<T,Index2D, SeparableRHS2D<T, Basis2D>, NoPrec2D>             AdaptRHS;
typedef AverageOutput2D<T, Index2D, Basis2D>							 Output;

// Algorithm Definition
    // Class for calling the truth solver for snapshot calculations
typedef CompressionWeightedPDE2D<T, Basis2D>                             		 Compression;
typedef IndexsetTruthSolver<T, Basis2D, DiagPrec2D, Index2D, Compression>   	 IndexsetSolver;
    // Class containing all \calN-dependent data and functions
typedef AdaptiveRBTruth2D<T, Basis2D, DiagPrec2D,IndexsetSolver, Compression>     RBTruth;
    // Class containing only N-dependent data and functions
typedef RBModel2D<T, RBTruth>                                           		  RBModel;

// Data type definitions
typedef Coefficients<Lexicographical,T,Index2D>                       Coeffs;
typedef flens::DenseVector<flens::Array<T> >                          DenseVectorT;
typedef flens::GeMatrix<flens::FullStorage<T, cxxblas::ColMajor> >    FullColMatrixT;


/* Example: Thermal block mit Testfeld
 *
 *          - theta * u_xx = f on (0,1)^2
 *         u(x,0) = u(x,1) = 0
 *
 *        Here: theta = theta_1 on (0, 0.5) x (0,1)
 *              theta = theta_2 on (0,5, 1) x (0,1)
 *                  f = 1       on (0.4, 0.6) x (0.45, 0.55)
 *                  omega_test = (0.7,0.7) x (0.8,0.8)
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
output_x(T x)
{
	return (x >= 0.7) && (x<= 0.8) ? 1 : 0 ;
}

T
output_y(T y)
{
	return (y >= 0.7) && (y<= 0.8) ? 1 : 0 ;
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
        cerr << "Usage " << argv[0] << " d d_ j0_x j0_y indexsetfile offline_data_directory outputfile_basename" << endl;
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
    bool use_A_operator_matrices = true;
    RBTruth rb_truth(basis2d, prec, use_inner_product_matrix, use_A_operator_matrices);
    rb_model.set_truthmodel(rb_truth);

        // Attach an inner product (here: H1 norm)
    AdaptHHOp2D h1norm(basis2d, 1., noprec, 1e-10);
    rb_truth.attach_inner_product_op(h1norm);



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

    WeightedAdaptHHOp2D hh_1(basis2d, 0, w_x_1, w_y, noprec, 1e-10);
    WeightedAdaptHHOp2D hh_2(basis2d, 0, w_x_2, w_y, noprec, 1e-10);

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


    /*RHS for Output Functional*/
    DenseVectorT singpts_output_x(4), singpts_output_y(4);
    singpts_output_x = 0., 0.7, 0.8, 1.;
    singpts_output_y = 0., 0.7, 0.8, 1.;
    SeparableFunction2D<T> sep_output_function(output_x, singpts_output_x, output_y, singpts_output_y);
    SeparableRHS2D<T, Basis2D> sep_output_rhs(basis2d, sep_output_function, noDeltas, noDeltas, 4);
    AdaptRHS Output_RHS(sep_output_rhs, noprec);
    Output AverageOutput(basis2d, 0.7,0.8,0.7,0.8, Output_RHS);

    rb_model.truth->attach_F_q(theta_f_1, F);

    cout << "Q_f = " << rb_model.Q_f() << endl;

    /* Truth Solver: on fixed indexset */
    IndexSet<Index2D> basisset;

    readIndexSet2D<T>(basisset, argv[5]);
    std::cout << "Basis Set: " << basisset.size() << " functions" << std::endl;

    SolverCall call = call_cg;
    IndexsetSolver indexset_solver(basisset, rb_truth, call);
    rb_truth.set_truthsolver(indexset_solver);
    rb_truth.assemble_inner_product_matrix(basisset);

    rb_model.read_RB_data(argv[6]);
    stringstream bf_offline_data;
    bf_offline_data << argv[6] << "/bf";
    rb_model.read_basis_functions(bf_offline_data.str().c_str());

    cout << "=================================================================" << endl;
    cout << "=====       ONLINE : TESTS                                 ======" << endl;
    cout << "=================================================================" << endl;

    string testerrorfile = argv[7];
    testerrorfile = testerrorfile + "_Test.txt";

    /* Test */

    // Construct test set
    std::vector<std::vector<T> > Xi_test;
    T test_min = 0.13;
    T test_max = 1.96;
    int n_test = 29;
    T h_test = (test_max - test_min) / (n_test-1);
    cout << "Testset: " << endl;
    for(int i = 0; i < n_test; ++i){
      vector<T> test_param;
      test_param.push_back(test_min + h_test * i);
      Xi_test.push_back(test_param);
      cout << Xi_test[i][0] << endl;
    }
    cout << endl;


    vector<T> maxerr(rb_model.n_bf());
    vector<T> maxerr_mu(rb_model.n_bf());
    vector<T> maxerr_bound(rb_model.n_bf());

    // Find maximum error over test set for all n = 1,...,N
    for(unsigned int i = 0; i < Xi_test.size(); ++i){
      cout << endl << "---- mu = " << Xi_test[i][0] << " --------------------------" << endl << endl;

      // Reference truth solve at test parameter
      rb_model.set_current_param(Xi_test[i]);
      Coeffs u = rb_model.truth->solver->truth_solve();
      cout << "size coeffs u : " << u.size() << endl;
      T output_u = AverageOutput.operator()(u);
      cout << "--------output operator u : " << output_u << endl;

      // RB solves for different basis sizes
      for(unsigned int n = 1; n <= rb_model.n_bf(); ++n){

        DenseVectorT u_N = rb_model.RB_solve(n, Xi_test[i]);
        T error_bound = rb_model.residual_dual_norm(u_N, Xi_test[i]) / rb_model.alpha_LB(Xi_test[i]);

        Coeffs u_approx = rb_model.reconstruct_u_N(u_N, n);
        cout << "size coeffs u_approx : " << u_approx.size() << endl<< endl;

        T output_u_approx = AverageOutput.operator()(u_approx);
        cout << "--------output operator u_approx : " << output_u_approx << endl;
        Coeffs coeff_diff = u - u_approx;
        T err_norm = rb_model.truth->inner_product(coeff_diff, coeff_diff);

        cout << "Differenz der Outputs: " << abs(output_u-output_u_approx) << endl;

        cout << "   N = " << n << ": " << err_norm << " " << error_bound << endl;

        if(err_norm > maxerr[n-1]){
          maxerr[n-1] = err_norm;
          maxerr_mu[n-1] = Xi_test[i][0];
          maxerr_bound[n-1] = error_bound;
        }
      }
    }





//    T output_u_diff = AverageOutput.operator()(coeff_diff);
  //  cout << "--------output operator u_diff : " << output_u_diff << endl;

    // Print to errorfile
    ofstream errorfile(testerrorfile.c_str());
    errorfile << "# N MaxError(H1Norm) ErrorBound Mu" << endl;
    for(unsigned int n = 0; n < rb_model.n_bf(); n++){
      errorfile << n+1 << " " << maxerr[n] << " " << maxerr_bound[n] << " " << maxerr_mu[n] << endl;
    }
    errorfile.close();

    return 0;
}
