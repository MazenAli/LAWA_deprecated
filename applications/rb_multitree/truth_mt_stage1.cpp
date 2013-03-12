#include <iostream>
#include <iomanip>
#include <utility>
#include <array>
#include <vector>
#include <lawa/lawa.h>

#include <lawa/methods/rb/datastructures/thetastructure.h>
#include <lawa/methods/rb/operators/affinelocaloperator.h>
#include <lawa/methods/rb/righthandsides/affinerhs.h>
#include <lawa/methods/adaptive/solvers/multitreeawgm2.h>
#include <lawa/methods/rb/righthandsides/flexiblebilformrhs.h>

using namespace std;
using namespace lawa;

//===============================================================//
//========= TYPEDEFS  =======================//
//===============================================================//

//==== General ====//
typedef double T;
typedef flens::GeMatrix<flens::FullStorage<T, cxxblas::ColMajor> >  FullColMatrixT;
typedef flens::DenseVector<flens::Array<T> >                        DenseVectorT;

//==== Basis 1D & 2D ====//
typedef Basis<T, Primal, Periodic, CDF>	                            TrialBasis_Time;
//typedef Basis<T, Primal, Periodic, CDF>								TestBasis_Time;
typedef Basis<T, Primal, Interval, Dijkema>						    TestBasis_Time;
typedef Basis<T, Primal, Interval, Dijkema>							Basis_Space;
typedef TrialBasis_Time::RefinementBasis                            TrialRefBasis_Time;
typedef TestBasis_Time::RefinementBasis                           	TestRefBasis_Time;
typedef Basis_Space::RefinementBasis                                RefBasis_Space;

// !!! ATTENTION: Use typedefs in space for definition of space-time inner product !!!
//    => Assumption that test basis in time = space basis

typedef TensorBasis2D<Adaptive,TrialBasis_Time,Basis_Space>         Basis2D_Trial;
typedef TensorBasis2D<Adaptive,TestBasis_Time,Basis_Space>          Basis2D_Test;

//==== Adaptive Operators ====//
typedef AdaptiveConvectionOperator1D_PG<T,TrialBasis_Time,TestBasis_Time>                   Convection1D_Time;
typedef AdaptiveIdentityOperator1D_PG<T,TrialBasis_Time,TestBasis_Time>                     Identity1D_Time;
typedef AdaptiveIdentityOperator1D_PG<T,Basis_Space,Basis_Space>                            Identity1D_Space;
typedef AdaptiveLaplaceOperator1D_PG<T,Basis_Space,Basis_Space>                             Laplace1D_Space;

typedef TransposedAdaptiveConvectionOperator1D_PG<T,TrialBasis_Time,TestBasis_Time>         TranspConvection1D_Time;
typedef TransposedAdaptiveIdentityOperator1D_PG<T,TrialBasis_Time,TestBasis_Time>           TranspIdentity1D_Time;
typedef TransposedAdaptiveIdentityOperator1D_PG<T,Basis_Space,Basis_Space>                  TranspIdentity1D_Space;
typedef TransposedAdaptiveLaplaceOperator1D_PG<T,Basis_Space,Basis_Space>                   TranspLaplace1D_Space;

typedef AdaptiveConvectionOperator1D_PG<T,TrialRefBasis_Time,TestRefBasis_Time>             RefConvection1D_Time;
typedef AdaptiveIdentityOperator1D_PG<T,TrialRefBasis_Time,TestRefBasis_Time>               RefIdentity1D_Time;
typedef AdaptiveIdentityOperator1D_PG<T,RefBasis_Space,RefBasis_Space>                      RefIdentity1D_Space;
typedef AdaptiveLaplaceOperator1D_PG<T,RefBasis_Space,RefBasis_Space>                       RefLaplace1D_Space;

typedef TransposedAdaptiveConvectionOperator1D_PG<T,TrialRefBasis_Time,TestRefBasis_Time>   RefTranspConvection1D_Time;
typedef TransposedAdaptiveIdentityOperator1D_PG<T,TrialRefBasis_Time,TestRefBasis_Time>     RefTranspIdentity1D_Time;
typedef TransposedAdaptiveIdentityOperator1D_PG<T,RefBasis_Space,RefBasis_Space>            RefTranspIdentity1D_Space;
typedef TransposedAdaptiveLaplaceOperator1D_PG<T,RefBasis_Space,RefBasis_Space>             RefTranspLaplace1D_Space;

//==== LocalOperators ====//
typedef LocalOperator1D<TestBasis_Time,TrialBasis_Time, 
                        RefConvection1D_Time,Convection1D_Time>                 LOp_Conv1D_Time;
typedef LocalOperator1D<TestBasis_Time,TrialBasis_Time, 
                        RefIdentity1D_Time,Identity1D_Time>                     LOp_Id1D_Time;                        
typedef LocalOperator1D<Basis_Space,Basis_Space,
						RefIdentity1D_Space,Identity1D_Space>	                LOp_Id1D_Space;
typedef LocalOperator1D<Basis_Space,Basis_Space,
						RefLaplace1D_Space,Laplace1D_Space>	                    LOp_Lapl1D_Space;
						
typedef LocalOperator1D<TrialBasis_Time,TestBasis_Time, 
                        RefTranspConvection1D_Time,TranspConvection1D_Time>     LOpT_Conv1D_Time;
typedef LocalOperator1D<TrialBasis_Time,TestBasis_Time, 
                        RefTranspIdentity1D_Time,TranspIdentity1D_Time>         LOpT_Id1D_Time;                        
typedef LocalOperator1D<Basis_Space,Basis_Space,
						RefTranspIdentity1D_Space,TranspIdentity1D_Space>	    LOpT_Id1D_Space;
typedef LocalOperator1D<Basis_Space,Basis_Space,
						RefTranspLaplace1D_Space,TranspLaplace1D_Space>	        LOpT_Lapl1D_Space;

typedef LocalOperator2D<LOp_Conv1D_Time, LOp_Id1D_Space>                        LOp_Conv_Id_2D;
typedef LocalOperator2D<LOp_Id1D_Time, LOp_Lapl1D_Space>                        LOp_Id_Lapl_2D;
typedef LocalOperator2D<LOpT_Conv1D_Time, LOpT_Id1D_Space>                      LOpT_Conv_Id_2D;
typedef LocalOperator2D<LOpT_Id1D_Time, LOpT_Lapl1D_Space>                      LOpT_Id_Lapl_2D;

typedef LocalOperator2D<LOp_Id1D_Space,LOp_Id1D_Space>							LOp_Test_Id_Id_2D;
typedef LocalOperator2D<LOp_Id1D_Space,LOp_Lapl1D_Space>						LOp_Test_Id_Lapl_2D;

//==== CompoundOperators ====//
typedef FlexibleCompoundLocalOperator<Index2D,AbstractLocalOperator2D<T> > 		Flex_COp_2D;
typedef CompoundLocalOperator<Index2D,LOp_Conv_Id_2D,LOp_Id_Lapl_2D>    		COp_Heat;
typedef CompoundLocalOperator<Index2D,LOpT_Conv_Id_2D,LOpT_Id_Lapl_2D>    	    COpT_Heat;
typedef CompoundLocalOperator<Index2D,LOp_Test_Id_Id_2D, LOp_Test_Id_Lapl_2D>	COp_TestInnProd;

//==== Preconditioners ====//
typedef AdaptiveLeftNormPreconditioner2D<T,Basis2D_Test>            LeftPrec2D;
typedef AdaptiveRightNormPreconditioner2D_c<T,Basis2D_Trial>        RightPrec2D;
typedef NoPreconditioner<T, Index2D>								NoPrec2D;

//==== RightHandSides ====//
typedef SeparableRHS2D<T,Basis2D_Test>                              SeparableRhsIntegral2D;
typedef RHS<T,Index2D,SeparableRhsIntegral2D,
            NoPrec2D>                                         		SeparableRhs;


//==== RB Stuff ====//
const size_t PDim = 1;
typedef AffineLocalOperator<Index2D,AbstractLocalOperator2D<T>,1>			Affine_Op_2D;
typedef AffineRhs<T,Index2D,SeparableRhs,1>									Affine_Rhs_2D;
typedef FlexibleCompoundRhs<T,Index2D,SeparableRhs>							RieszF_Rhs_2D;
typedef FlexibleBilformRhs<Index2D,AbstractLocalOperator2D<T> >				RieszA_Rhs_2D;

typedef MultiTreeAWGM_PG<Index2D,Basis2D_Trial, Basis2D_Test,Affine_Op_2D,
		Affine_Op_2D,Affine_Rhs_2D,RightPrec2D,LeftPrec2D>					MT_AWGM_Truth;
typedef MultiTreeAWGM2<Index2D,Basis2D_Test, COp_TestInnProd,
		RieszF_Rhs_2D,LeftPrec2D>											MT_AWGM_Riesz_F;
typedef MultiTreeAWGM2<Index2D,Basis2D_Test, COp_TestInnProd,
		RieszA_Rhs_2D,LeftPrec2D>											MT_AWGM_Riesz_A;


T f_t(T t)
{
    return std::cos(2*M_PI*t);
}

T f_x(T x)
{
    return -4*(x-0.5)*(x-0.5) + 1;
}

T dummy(T, T)
{
    return 0;
}

T no_theta(const std::array<T,PDim>& /*mu*/)
{
	return 1.;
}

T theta_diff(const std::array<T,PDim>& mu)
{
	return mu[0];
}


int main () {

	//===============================================================//
	//========= PROBLEM SETUP  =======================//
	//===============================================================//
	
    int d   = 2;
    int d_  = 2;
    int j0  = 2;

    //getchar();

    /// Basis initialization
    TrialBasis_Time      basis_per(d,d_,j0);
    TestBasis_Time       basis_int(d,d_,j0);
    Basis_Space 		 basis_intbc(d,d_,j0);
    basis_intbc.enforceBoundaryCondition<DirichletBC>();

    Basis2D_Trial basis2d_trial(basis_per,basis_intbc);
    Basis2D_Test  basis2d_test(basis_int,basis_intbc);

    /// Initialization of operator

    // Bilinear Forms
    Convection1D_Time			ConvectionBil_t(basis_per, basis_int);
    Identity1D_Time 		    IdentityBil_t(basis_per, basis_int);
    Identity1D_Space 	        IdentityBil_x(basis_intbc, basis_intbc);
    Laplace1D_Space 	        LaplaceBil_x(basis_intbc, basis_intbc);
    Identity1D_Space 	        IdentityBil_Test_t(basis_int, basis_int);
    
    RefConvection1D_Time 		RefConvectionBil_t(basis_per.refinementbasis, basis_int.refinementbasis);
    RefIdentity1D_Time 		    RefIdentityBil_t(basis_per.refinementbasis, basis_int.refinementbasis);
    RefIdentity1D_Space 	    RefIdentityBil_x(basis_intbc.refinementbasis, basis_intbc.refinementbasis);
    RefLaplace1D_Space 	        RefLaplaceBil_x(basis_intbc.refinementbasis, basis_intbc.refinementbasis);
    RefIdentity1D_Space 	    RefIdentityBil_Test_t(basis_int.refinementbasis, basis_int.refinementbasis);

    // Transposed Bilinear Forms
    TranspConvection1D_Time 	TranspConvectionBil_t(basis_per, basis_int);
    TranspIdentity1D_Time 		TranspIdentityBil_t(basis_per, basis_int);
    TranspIdentity1D_Space 	    TranspIdentityBil_x(basis_intbc, basis_intbc);
    TranspLaplace1D_Space 	    TranspLaplaceBil_x(basis_intbc, basis_intbc);
    
    RefTranspConvection1D_Time 	RefTranspConvectionBil_t(basis_per.refinementbasis, basis_int.refinementbasis);
    RefTranspIdentity1D_Time 	RefTranspIdentityBil_t(basis_per.refinementbasis, basis_int.refinementbasis);
    RefTranspIdentity1D_Space 	RefTranspIdentityBil_x(basis_intbc.refinementbasis, basis_intbc.refinementbasis);
    RefTranspLaplace1D_Space 	RefTranspLaplaceBil_x(basis_intbc.refinementbasis, basis_intbc.refinementbasis);

    /// Initialization of local operator
    LOp_Conv1D_Time		lOp_Conv1D_t	(basis_int, basis_per, RefConvectionBil_t, ConvectionBil_t);
    LOp_Id1D_Time		lOp_Id1D_t  	(basis_int, basis_per, RefIdentityBil_t, IdentityBil_t);
    LOp_Id1D_Space		lOp_Id1D_x  	(basis_intbc, basis_intbc, RefIdentityBil_x, IdentityBil_x);
    LOp_Lapl1D_Space	lOp_Lapl1D_x	(basis_intbc, basis_intbc, RefLaplaceBil_x, LaplaceBil_x);
    LOp_Id1D_Space		lOp_Id1D_Test_t (basis_int, basis_int, RefIdentityBil_Test_t, IdentityBil_Test_t);
    
    LOpT_Conv1D_Time	lOpT_Conv1D_t(basis_per, basis_int, RefTranspConvectionBil_t, TranspConvectionBil_t);
    LOpT_Id1D_Time		lOpT_Id1D_t  (basis_per, basis_int, RefTranspIdentityBil_t, TranspIdentityBil_t);
    LOpT_Id1D_Space		lOpT_Id1D_x  (basis_intbc, basis_intbc, RefTranspIdentityBil_x, TranspIdentityBil_x);
    LOpT_Lapl1D_Space	lOpT_Lapl1D_x(basis_intbc, basis_intbc, RefTranspLaplaceBil_x, TranspLaplaceBil_x);

    LOp_Conv_Id_2D		localConvectionIdentityOp2D		(lOp_Conv1D_t, 		lOp_Id1D_x);
    LOp_Id_Lapl_2D		localIdentityLaplaceOp2D		(lOp_Id1D_t, 		lOp_Lapl1D_x);
    LOp_Test_Id_Id_2D	localIdentityIdentityOp2D_Test	(lOp_Id1D_Test_t, 	lOp_Id1D_x);
    LOp_Test_Id_Lapl_2D	localIdentityLaplaceOp2D_Test	(lOp_Id1D_Test_t, 	lOp_Lapl1D_x);

    LOpT_Conv_Id_2D		transpLocalConvectionIdentityOp2D	(lOpT_Conv1D_t, lOpT_Id1D_x);
    LOpT_Id_Lapl_2D		transpLocalIdentityLaplaceOp2D		(lOpT_Id1D_t, lOpT_Lapl1D_x);

    localConvectionIdentityOp2D.setJ(9);
    localIdentityLaplaceOp2D.setJ(9);
    transpLocalConvectionIdentityOp2D.setJ(9);
    transpLocalIdentityLaplaceOp2D.setJ(9);
    localIdentityIdentityOp2D_Test.setJ(9);
    localIdentityLaplaceOp2D_Test.setJ(9);


    /// Initialization of preconditioner
    LeftPrec2D leftPrec(basis2d_test);
    RightPrec2D rightPrec(basis2d_trial);

    NoPrec2D noPrec;

    /// Initialization of rhs

    /// Right Hand Side:
    ///     No Singular Supports in both dimensions
    DenseVectorT sing_support;
    ///      Forcing Functions
    SeparableFunction2D<T> F1(f_t, sing_support, f_x, sing_support);
    ///     Peaks: points and corresponding coefficients
    ///             (heights of jumps in derivatives)
    FullColMatrixT nodeltas;
    SeparableRhsIntegral2D			rhs(basis2d_test, F1, nodeltas, nodeltas, 20);
    SeparableRhs           			F(rhs,noPrec);


	//===============================================================//
	//===============  RB SETUP =====================================//
	//===============================================================//

    // Affine Decompositions:
    // 	Left Hand Side
    vector<ThetaStructure<T,PDim>::ThetaFct>	lhs_theta_fcts;
    lhs_theta_fcts.push_back(no_theta);
    lhs_theta_fcts.push_back(theta_diff);
    ThetaStructure<T,PDim> lhs_theta(lhs_theta_fcts);

    vector<AbstractLocalOperator2D<T>* > lhs_ops, lhs_opsT;
    lhs_ops.push_back(&localConvectionIdentityOp2D);
    lhs_ops.push_back(&localIdentityLaplaceOp2D);
    lhs_opsT.push_back(&transpLocalConvectionIdentityOp2D);
    lhs_opsT.push_back(&transpLocalIdentityLaplaceOp2D);

    Affine_Op_2D affine_lhs(lhs_theta, lhs_ops);
    Affine_Op_2D affine_lhs_T(lhs_theta, lhs_opsT);

    COp_TestInnProd innprod_Y(localIdentityIdentityOp2D_Test, localIdentityLaplaceOp2D_Test);

    // Right Hand Side
    vector<ThetaStructure<T,PDim>::ThetaFct> rhs_theta_fcts;
    rhs_theta_fcts.push_back(no_theta);
    ThetaStructure<T,PDim> rhs_theta(rhs_theta_fcts);
    vector<SeparableRhs*> rhs_fcts;
    rhs_fcts.push_back(&F);

    Affine_Rhs_2D affine_rhs(rhs_theta, rhs_fcts);
    RieszF_Rhs_2D rieszF_rhs(rhs_fcts);
    RieszA_Rhs_2D rieszA_rhs(lhs_ops);

    std::array<T,PDim> mu = {{1.}};
    lhs_theta.set_param(mu);
    rhs_theta.set_param(mu);


	//===============================================================//
	//===============  AWGM =========================================//
	//===============================================================//


    //----------- Truth Solver ---------------- //

    /* AWGM PG Parameters Default Values
    double tol = 5e-03;
	double alpha = 0.7;
	size_t max_its = 100;
	size_t max_basissize = 400000;
	bool reset_res = false;
	bool print_info = true;
	bool verbose = true;
	bool plot_solution = false;
	bool verbose_extra = false; //(print added wavelet indizes)
	size_t hashmapsize_trial = 10;
	size_t hashmapsize_test = 10;
	std::string info_filename = "awgm_cgls_conv_info.txt";
	std::string plot_filename = "awgm_cgls_u_plot";
	*/

    /* CGLS Parameters Default Values
	bool adaptive_tol = true;
	size_t max_its = 100;
	double init_tol = 0.001;
	double res_reduction = 0.01;
	double absolute_tol = 1e-8;
	bool verbose = true;
	*/

    cout << "||=====================================================================||" << endl;
    cout << "||================  TRUTH SOLUTION ====================================||" << std::endl;
    cout << "||=====================================================================||" << endl << endl;

    AWGM_PG_Parameters awgm_truth_parameters;
    IS_Parameters cgls_parameters;
    awgm_truth_parameters.plot_solution = true;

    MT_AWGM_Truth multitree_awgm(basis2d_trial, basis2d_test, affine_lhs, affine_lhs_T,
    							 affine_rhs, rightPrec, leftPrec, awgm_truth_parameters, cgls_parameters);

    multitree_awgm.awgm_params.print();
    multitree_awgm.is_params.print();
    multitree_awgm.set_sol(dummy);

    /// Initialization of solution vector and initial index sets
    Coefficients<Lexicographical,T,Index2D> u;

    T gamma = 0.2;
    IndexSet<Index2D> LambdaTrial, LambdaTest;
    getSparseGridIndexSet(basis2d_trial,LambdaTrial,2,0,gamma);
    getSparseGridIndexSet(basis2d_test ,LambdaTest ,2,1,gamma);

    LambdaTrial.clear();
    LambdaTest.clear();
    getSparseGridIndexSet(basis2d_trial,LambdaTrial,2,0,gamma);
    getSparseGridIndexSet(basis2d_test ,LambdaTest ,2,1,gamma);

    Timer time;
    time.start();
    multitree_awgm.solve(u, LambdaTrial, LambdaTest);
    time.stop();
    cout << "Solution took " << time.elapsed() << " seconds" << endl;

    //----------- RieszF Solver ---------------- //

    /* AWGM Parameters Default Values
    double tol = 5e-03;
	double alpha = 0.7;
	size_t max_its = 100;
	size_t max_basissize = 400000;
	bool print_info = true;
	bool verbose = true;
	bool plot_solution = false;
	bool verbose_extra = false; //(print added wavelet indizes)
	size_t hashmapsize = 10;
	std::string info_filename = "awgm_cg_conv_info.txt";
	std::string plot_filename = "awgm_cg_u_plot";
	*/
    cout << "||=====================================================================||" << endl;
    cout << "||================  RIESZ SOLUTION F ==================================||" << std::endl;
    cout << "||=====================================================================||" << endl << endl;

    AWGM_Parameters awgm_riesz_f_parameters;
    awgm_riesz_f_parameters.plot_solution = true;
    awgm_riesz_f_parameters.info_filename = "awgm_R_f_conv_info.txt";
    awgm_riesz_f_parameters.plot_filename = "awgm_R_f_plot";

    MT_AWGM_Riesz_F multitree_awgm_rieszF(basis2d_test, innprod_Y, rieszF_rhs, leftPrec, awgm_riesz_f_parameters, cgls_parameters);

    multitree_awgm_rieszF.awgm_params.print();
    multitree_awgm_rieszF.is_params.print();
    multitree_awgm_rieszF.set_sol(dummy);

    Coefficients<Lexicographical,T,Index2D> r_f;

    LambdaTest.clear();
    getSparseGridIndexSet(basis2d_test ,LambdaTest ,2,1,gamma);

    time.start();
    //multitree_awgm_rieszF.solve(r_f, LambdaTest);
    time.stop();
    cout << "Solution took " << time.elapsed() << " seconds" << endl;

    //----------- RieszA Solver ---------------- //

    cout << "||=====================================================================||" << endl;
    cout << "||================  RIESZ SOLUTION A ==================================||" << std::endl;
    cout << "||=====================================================================||" << endl << endl;

    AWGM_Parameters awgm_riesz_a_parameters;
    awgm_riesz_a_parameters.tol = 5e-04;
    awgm_riesz_a_parameters.plot_solution = true;

    MT_AWGM_Riesz_A multitree_awgm_rieszA(basis2d_test, innprod_Y, rieszA_rhs, leftPrec, awgm_riesz_a_parameters, cgls_parameters);

    multitree_awgm_rieszA.set_sol(dummy);

    Coefficients<Lexicographical,T,Index2D> r_a_0, r_a_1;

    LambdaTest.clear();
    getSparseGridIndexSet(basis2d_test ,LambdaTest ,2,1,gamma);

    rieszA_rhs.set_active_u(&u);

    rieszA_rhs.set_active_comp(0);
    multitree_awgm_rieszA.awgm_params.info_filename = "awgm_R_a_0_conv_info.txt";
    multitree_awgm_rieszA.awgm_params.plot_filename = "awgm_R_a_0_plot";
    multitree_awgm_rieszA.awgm_params.print();
    multitree_awgm_rieszA.is_params.print();
    time.start();
    //multitree_awgm_rieszA.solve(r_a_0, LambdaTest);
    time.stop();
    cout << "Solution took " << time.elapsed() << " seconds" << endl;

    cout << "||=====================================================================||" << endl << endl;

    LambdaTest.clear();
    getSparseGridIndexSet(basis2d_test ,LambdaTest ,2,1,gamma);

    rieszA_rhs.set_active_comp(1);
    multitree_awgm_rieszA.awgm_params.info_filename = "awgm_R_a_1_conv_info.txt";
    multitree_awgm_rieszA.awgm_params.plot_filename = "awgm_R_a_1_plot";
    multitree_awgm_rieszA.awgm_params.print();
    multitree_awgm_rieszA.is_params.print();
    time.start();
    multitree_awgm_rieszA.solve(r_a_1, LambdaTest);
    time.stop();
    cout << "Solution took " << time.elapsed() << " seconds" << endl;

    return 0;
}

