#include <iostream>
#include <iomanip>
#include <utility>
#include <array>
#include <lawa/lawa.h>

#include <lawa/methods/rb/datastructures/thetastructure.h>
#include <lawa/methods/rb/operators/affinelocaloperator.h>
#include <lawa/methods/rb/righthandsides/affinerhs.h>


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

//==== CompoundOperators ====//
typedef FlexibleCompoundLocalOperator<Index2D,AbstractLocalOperator2D<T> > 		Flex_COp_2D;
typedef CompoundLocalOperator<Index2D,LOp_Conv_Id_2D,LOp_Id_Lapl_2D>    		COp_Heat;
typedef CompoundLocalOperator<Index2D,LOpT_Conv_Id_2D,LOpT_Id_Lapl_2D>    	    COpT_Heat;

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
typedef MultiTreeAWGM_PG<Index2D,Basis2D_Trial, Basis2D_Test,Affine_Op_2D,
		Affine_Op_2D,Affine_Rhs_2D,RightPrec2D,LeftPrec2D>				MT_AWGM;

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

T no_theta(const std::array<T,PDim>& mu)
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
    
    RefConvection1D_Time 		RefConvectionBil_t(basis_per.refinementbasis, basis_int.refinementbasis);
    RefIdentity1D_Time 		    RefIdentityBil_t(basis_per.refinementbasis, basis_int.refinementbasis);
    RefIdentity1D_Space 	    RefIdentityBil_x(basis_intbc.refinementbasis, basis_intbc.refinementbasis);
    RefLaplace1D_Space 	        RefLaplaceBil_x(basis_intbc.refinementbasis, basis_intbc.refinementbasis);

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
    LOp_Conv1D_Time		lOp_Conv1D_t(basis_int, basis_per, RefConvectionBil_t, ConvectionBil_t);
    LOp_Id1D_Time		lOp_Id1D_t  (basis_int, basis_per, RefIdentityBil_t, IdentityBil_t);
    LOp_Id1D_Space		lOp_Id1D_x  (basis_intbc, basis_intbc, RefIdentityBil_x, IdentityBil_x);
    LOp_Lapl1D_Space	lOp_Lapl1D_x(basis_intbc, basis_intbc, RefLaplaceBil_x, LaplaceBil_x);
    
    LOpT_Conv1D_Time	lOpT_Conv1D_t(basis_per, basis_int, RefTranspConvectionBil_t, TranspConvectionBil_t);
    LOpT_Id1D_Time		lOpT_Id1D_t  (basis_per, basis_int, RefTranspIdentityBil_t, TranspIdentityBil_t);
    LOpT_Id1D_Space		lOpT_Id1D_x  (basis_intbc, basis_intbc, RefTranspIdentityBil_x, TranspIdentityBil_x);
    LOpT_Lapl1D_Space	lOpT_Lapl1D_x(basis_intbc, basis_intbc, RefTranspLaplaceBil_x, TranspLaplaceBil_x);

    LOp_Conv_Id_2D		localConvectionIdentityOp2D(lOp_Conv1D_t, lOp_Id1D_x);
    LOp_Id_Lapl_2D		localIdentityLaplaceOp2D(lOp_Id1D_t, lOp_Lapl1D_x);
    
    LOpT_Conv_Id_2D		transpLocalConvectionIdentityOp2D(lOpT_Conv1D_t, lOpT_Id1D_x);
    LOpT_Id_Lapl_2D		transpLocalIdentityLaplaceOp2D(lOpT_Id1D_t, lOpT_Lapl1D_x);

    localConvectionIdentityOp2D.setJ(9);
    localIdentityLaplaceOp2D.setJ(9);
    transpLocalConvectionIdentityOp2D.setJ(9);
    transpLocalIdentityLaplaceOp2D.setJ(9);

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
    vector<ThetaStructure<T,PDim>::ThetaFct>	lhs_thetas;
    lhs_thetas.push_back(no_theta);
    lhs_thetas.push_back(theta_diff);
    vector<AbstractLocalOperator2D<T>* > lhs_ops, lhs_opsT;
    lhs_ops.push_back(&localConvectionIdentityOp2D);
    lhs_ops.push_back(&localIdentityLaplaceOp2D);
    lhs_opsT.push_back(&transpLocalConvectionIdentityOp2D);
    lhs_opsT.push_back(&transpLocalIdentityLaplaceOp2D);

    Affine_Op_2D affine_lhs(lhs_thetas, lhs_ops);
    Affine_Op_2D affine_lhs_T(lhs_thetas, lhs_opsT);

    // Right Hand Side
    vector<ThetaStructure<T,PDim>::ThetaFct> rhs_thetas;
    rhs_thetas.push_back(no_theta);
    vector<SeparableRhs*> rhs_fcts;
    rhs_fcts.push_back(&F);

    Affine_Rhs_2D affine_rhs(rhs_thetas, rhs_fcts);

	//===============================================================//
	//===============  AWGM =========================================//
	//===============================================================//


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
	*/

    /* CGLS Parameters Default Values
	bool adaptive_tol = true;
	size_t max_its = 100;
	double init_tol = 0.001;
	double res_reduction = 0.01;
	double absolute_tol = 1e-8;
	bool verbose = true;
	*/

    // If you want other parameters
    AWGM_PG_Parameters awgm_parameters;
    IS_Parameters cgls_parameters;
    // .... set them here:
    awgm_parameters.plot_solution = true;

    MT_AWGM multitree_awgm(basis2d_trial, basis2d_test, affine_lhs, affine_lhs_T,
    						affine_rhs, rightPrec, leftPrec, awgm_parameters, cgls_parameters);


    multitree_awgm.awgm_params.print();
    multitree_awgm.is_params.print();

    multitree_awgm.set_sol(dummy);

    /// Initialization of solution vector and initial index sets
    Coefficients<Lexicographical,T,Index2D> u;

    T gamma = 0.2;
    IndexSet<Index2D> LambdaTrial, LambdaTest;
    getSparseGridIndexSet(basis2d_trial,LambdaTrial,2,0,gamma);
    getSparseGridIndexSet(basis2d_test ,LambdaTest ,2,1,gamma);

    Timer time;
    time.start();
    multitree_awgm.solve(u, LambdaTrial, LambdaTest);
    time.stop();
    cout << "Solution took " << time.elapsed() << " seconds" << endl;

    return 0;
}

