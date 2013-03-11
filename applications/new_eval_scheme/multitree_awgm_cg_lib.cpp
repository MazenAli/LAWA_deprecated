#include <iostream>
#include <iomanip>
#include <utility>
#include <lawa/lawa.h>
#include <lawa/methods/adaptive/operators/localoperators/abstractlocaloperator2d.h>
#include <lawa/methods/adaptive/operators/localoperators/flexiblecompoundlocaloperator.h>
#include <lawa/methods/adaptive/solvers/multitreeawgm_pg.h>
#include <lawa/methods/adaptive/solvers/multitreeawgm2.h>
#include <lawa/methods/adaptive/algorithms/indexset_generation.h>
#include <lawa/methods/adaptive/solvers/solver_parameters.h>


using namespace std;
using namespace lawa;

typedef double T;

typedef flens::GeMatrix<flens::FullStorage<T, cxxblas::ColMajor> >   FullColMatrixT;
typedef flens::DenseVector<flens::Array<T> >                         DenseVectorT;

typedef Basis<T, Primal, Periodic, CDF>	                            TrialBasis_Time;
//typedef Basis<T, Primal, Periodic, CDF>								TestBasis_Time;
typedef Basis<T, Primal, Interval, Dijkema>							TestBasis_Time;
typedef Basis<T, Primal, Interval, Dijkema>							Basis_Space;
typedef TrialBasis_Time::RefinementBasis                            TrialBasis_Time_RefinementBasis;
typedef TestBasis_Time::RefinementBasis                           	TestBasis_Time_RefinementBasis;
typedef Basis_Space::RefinementBasis                                Basis_Space_RefinementBasis;

typedef TensorBasis2D<Adaptive,TrialBasis_Time,Basis_Space>             Basis2D_Trial;
typedef TensorBasis2D<Adaptive,TestBasis_Time,Basis_Space>             	Basis2D_Test;

typedef AdaptiveWeightedPDEOperator1D_PG<T,TrialBasis_Time,TestBasis_Time>    			PDEBilinearForm1D_Time;
typedef AdaptiveWeightedPDEOperator1D_PG<T,Basis_Space,Basis_Space>   					PDEBilinearForm1D_Space;
typedef TransposedAdaptiveWeightedPDEOperator1D_PG<T,TrialBasis_Time,TestBasis_Time>    TransPDEBilinearForm1D_Time;
typedef TransposedAdaptiveWeightedPDEOperator1D_PG<T,Basis_Space,Basis_Space>   		TransPDEBilinearForm1D_Space;

typedef AdaptiveWeightedPDEOperator1D_PG<T,
			TrialBasis_Time_RefinementBasis,TestBasis_Time_RefinementBasis>    			PDERefinementBilinearForm1D_Time;
typedef AdaptiveWeightedPDEOperator1D_PG<T,
			Basis_Space_RefinementBasis,Basis_Space_RefinementBasis>    				PDERefinementBilinearForm1D_Space;

typedef TransposedAdaptiveWeightedPDEOperator1D_PG<T,
			TrialBasis_Time_RefinementBasis,TestBasis_Time_RefinementBasis>    			TransPDERefinementBilinearForm1D_Time;
typedef TransposedAdaptiveWeightedPDEOperator1D_PG<T,
			Basis_Space_RefinementBasis,Basis_Space_RefinementBasis>    				TransPDERefinementBilinearForm1D_Space;

typedef LocalOperator1D<TestBasis_Time,TrialBasis_Time,
						PDERefinementBilinearForm1D_Time,
						PDEBilinearForm1D_Time>			        	LocalOp1D_t;
typedef LocalOperator1D<Basis_Space,Basis_Space,
						PDERefinementBilinearForm1D_Space,
						PDEBilinearForm1D_Space>			        LocalOp1D_x;
typedef LocalOperator1D<TrialBasis_Time,TestBasis_Time,
						TransPDERefinementBilinearForm1D_Time,
						TransPDEBilinearForm1D_Time>			    TransLocalOp1D_t;
typedef LocalOperator1D<Basis_Space,Basis_Space,
						TransPDERefinementBilinearForm1D_Space,
						TransPDEBilinearForm1D_Space>			    TransLocalOp1D_x;
typedef LocalOperator2D<LocalOp1D_t, LocalOp1D_x>                 	LocalOp2D;
typedef LocalOperator2D<TransLocalOp1D_t, TransLocalOp1D_x>         TransLocalOp2D;

typedef LocalOperator2D<LocalOp1D_x,LocalOp1D_x>					LocalOp2D_Test;

typedef FlexibleCompoundLocalOperator<Index2D,AbstractLocalOperator2D<T> > 	FlexibleCompoundLocalOperator2D;
typedef CompoundLocalOperator<Index2D,LocalOp2D,LocalOp2D>    				CompoundLocalOperator2D;
typedef CompoundLocalOperator<Index2D,TransLocalOp2D,TransLocalOp2D>    	TransCompoundLocalOperator2D;
typedef CompoundLocalOperator<Index2D,LocalOp2D_Test,LocalOp2D_Test>		CompoundLocalOperator2D_Test;


typedef AdaptiveLeftNormPreconditioner2D<T,Basis2D_Test>            LeftPrec2D;
typedef AdaptiveRightNormPreconditioner2D_c<T,Basis2D_Trial>        RightPrec2D;

typedef NoPreconditioner<T, Index2D>								NoPrec2D;

//Righthandsides definitions (separable)
typedef SeparableRHS2D<T,Basis2D_Test>                              SeparableRhsIntegral2D;

typedef SumOfTwoRHSIntegrals<T,Index2D,SeparableRhsIntegral2D,
                             SeparableRhsIntegral2D>                SumOfSeparableRhsIntegral2D;
typedef RHS<T,Index2D,SumOfSeparableRhsIntegral2D,
            NoPrec2D>                                         		SumOfSeparableRhs;

typedef MultiTreeAWGM_PG<Index2D,Basis2D_Trial, Basis2D_Test,CompoundLocalOperator2D,
			TransCompoundLocalOperator2D,SumOfSeparableRhs,RightPrec2D,LeftPrec2D>				MT_AWGM;
//typedef MultiTreeAWGM_PG<Index2D,Basis2D_Trial, Basis2D_Test,FlexibleCompoundLocalOperator2D,
//		FlexibleCompoundLocalOperator2D,SumOfSeparableRhs,RightPrec2D,LeftPrec2D>				MT_AWGM;

typedef MultiTreeAWGM2<Index2D,Basis2D_Test,CompoundLocalOperator2D_Test,
						SumOfSeparableRhs,LeftPrec2D>											MT_AWGM_Galerkin;

void
readIndexSetFromFile(IndexSet<Index2D> &Lambda, string filename);


T u1(T t)
{
    return std::cos(2*M_PI*t);
}

T u2(T x)
{
    return -4*(x-0.5)*(x-0.5) + 1;
}

T f_rhs_t(T t)
{
    return -2*M_PI* std::sin(2*M_PI*t);
}

T f_rhs_x(T /*y*/)
{
    return 8.;
}

T dummy(T ,T)
{
    return 0;
}

T zero_fct(T /*x*/){
	return 0;
}

T one_fct(T /*x*/){
	return 1;
}


int main () {

	//===============================================================//
	//========= PROBLEM SETUP  =======================//
	//===============================================================//

    /*if (argc!=5) {
        cout << "Usage: " << argv[0] << " d d_ j0 J" << endl;
        return 0;
    }*/
    int d   = 2;
    int d_  = 2;
    int j0  = 2;

    //getchar();

    /// Basis initialization
    TestBasis_Time       basis_int(d,d_,j0);
    Basis_Space 		 basis_intbc(d,d_,j0);
    basis_intbc.enforceBoundaryCondition<DirichletBC>();

    Basis2D_Test  basis2d_test(basis_int,basis_intbc);

    /// Initialization of operator
//    AdaptiveWeightedPDEOperator1D(const Basis1D& _basis1d, Function<T> &_reaction_f,
//                                  Function<T> &_convection_f, Function<T>& _diffusion_f,
//                                  int order=10,
//                                  bool reactionIsZero=false, bool convectionIsZero=false,
//                                  bool diffusionIsZero=false);
    DenseVectorT no_singPts;
    Function<T> zero_Fct(zero_fct,no_singPts);
    Function<T> one_Fct(one_fct,no_singPts);

    // Bilinear Forms
    PDEBilinearForm1D_Space 	IdentityBil_x(basis_intbc, basis_intbc, one_Fct, zero_Fct, zero_Fct, 10, false, true, true);
    PDEBilinearForm1D_Space 	LaplaceBil_x(basis_intbc, basis_intbc, zero_Fct, zero_Fct, one_Fct, 10, true, true, false);
    PDEBilinearForm1D_Space 	IdentityBil_Test_t(basis_int, basis_int, one_Fct, zero_Fct, zero_Fct, 10, false, true, true);


    PDERefinementBilinearForm1D_Space 	RefIdentityBil_x(basis_intbc.refinementbasis, basis_intbc.refinementbasis, one_Fct, zero_Fct, zero_Fct, 10, false, true, true);
    PDERefinementBilinearForm1D_Space 	RefLaplaceBil_x(basis_intbc.refinementbasis, basis_intbc.refinementbasis, zero_Fct, zero_Fct, one_Fct, 10, true, true, false);
    PDERefinementBilinearForm1D_Space 	RefIdentityBil_Test_t(basis_int.refinementbasis, basis_int.refinementbasis, one_Fct, zero_Fct, zero_Fct, 10, false, true, true);


    /// Initialization of local operator
    LocalOp1D_x		localIdentityOp1D_x(basis_intbc, basis_intbc, RefIdentityBil_x, IdentityBil_x);
    LocalOp1D_x		localLaplaceOp1D_x(basis_intbc, basis_intbc, RefLaplaceBil_x, LaplaceBil_x);
    LocalOp1D_x		localIdentityOp1D_Test_t(basis_int, basis_int, RefIdentityBil_Test_t, IdentityBil_Test_t);

    LocalOp2D_Test 		localIdentityIdentityOp2D_Test(localIdentityOp1D_Test_t, localIdentityOp1D_x);
    LocalOp2D_Test 		localIdentityLaplaceOp2D_Test(localIdentityOp1D_Test_t, localLaplaceOp1D_x);

    localIdentityIdentityOp2D_Test.setJ(9);
    localIdentityLaplaceOp2D_Test.setJ(9);

    // Use CompoundLocalOperator2D
    CompoundLocalOperator2D_Test       localOperator2D(localIdentityIdentityOp2D_Test,localIdentityLaplaceOp2D_Test);

    // Use FlexibleCompoundLocalOperator2D
//    vector<AbstractLocalOperator2D<T>* > localOperatorVec, transLocalOperatorVec;
//    localOperatorVec.push_back(&localConvectionIdentityOp2D);
//    localOperatorVec.push_back(&localIdentityLaplaceOp2D);
//    transLocalOperatorVec.push_back(&transLocalConvectionIdentityOp2D);
//    transLocalOperatorVec.push_back(&transLocalIdentityLaplaceOp2D);
//    FlexibleCompoundLocalOperator2D       localOperator2D(localOperatorVec);
//    FlexibleCompoundLocalOperator2D  	  transLocalOperator2D(transLocalOperatorVec);

    /// Initialization of preconditioner
    LeftPrec2D leftPrec(basis2d_test);

    NoPrec2D noPrec;

    /// Initialization of rhs

    /// Right Hand Side:
    ///     No Singular Supports in both dimensions
    DenseVectorT sing_support;
    ///      Forcing Functions
    SeparableFunction2D<T> F1(f_rhs_t, sing_support, u2, sing_support);
    SeparableFunction2D<T> F2(u1, sing_support, f_rhs_x, sing_support);
    ///     Peaks: points and corresponding coefficients
    ///             (heights of jumps in derivatives)
    FullColMatrixT nodeltas;
    SeparableRhsIntegral2D			rhs1(basis2d_test, F1, nodeltas, nodeltas, 20);
    SeparableRhsIntegral2D 			rhs2(basis2d_test, F2, nodeltas, nodeltas, 20);
    SumOfSeparableRhsIntegral2D 	rhsintegral2d(rhs1, rhs2);
    SumOfSeparableRhs           	F(rhsintegral2d,noPrec);


	//===============================================================//
	//===============  AWGM2 =========================================//
	//===============================================================//

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
	*/

    /* CG Parameters Default Values
	bool adaptive_tol = true;
	size_t max_its = 100;
	double init_tol = 0.001;
	double res_reduction = 0.01;
	double absolute_tol = 1e-8;
	bool verbose = true;
	*/

    AWGM_Parameters awgm_parameters;
    IS_Parameters	cg_parameters;
    awgm_parameters.plot_solution = true;

    MT_AWGM_Galerkin multitree_awgm_cg(basis2d_test, localOperator2D, F, leftPrec, awgm_parameters, cg_parameters);


    multitree_awgm_cg.awgm_params.print();
    multitree_awgm_cg.is_params.print();

    multitree_awgm_cg.set_sol(dummy);

    /// Initialization of solution vector and initial index sets
    Coefficients<Lexicographical,T,Index2D> u;

    T gamma = 0.2;
    IndexSet<Index2D> LambdaTest;
    getSparseGridIndexSet(basis2d_test ,LambdaTest ,2,1,gamma);

    Timer time;
    time.start();
    multitree_awgm_cg.solve(u, LambdaTest);
    time.stop();
    cout << "Solution took " << time.elapsed() << " seconds" << endl;

    return 0;
}

void
readIndexSetFromFile(IndexSet<Index2D> &Lambda, string filename)
{
    std::ifstream infile (filename.c_str());
    if (infile.is_open()) {
        cerr << "   Indexset file is open." << endl;
    }
    else {
        cerr << "   Indexset file " << filename.c_str()  << " is not open." << endl;
    }

	int t1,t2;
    int j1,j2;
    long k1,k2;
    T coeff;

    while(!infile.eof()) {

    	infile >> t1 >> j1 >> k1 >> t2 >> j2 >> k2 >> coeff;

        if (t1 == 1 && t2 == 1) {
            Index1D index_x(j1,k1,XWavelet);
            Index1D index_y(j2,k2,XWavelet);
            Lambda.insert(Index2D(index_x,index_y));
        }
        else if (t1 == 1 && t2 == 0) {
            Index1D index_x(j1,k1,XWavelet);
            Index1D index_y(j2,k2,XBSpline);
            Lambda.insert(Index2D(index_x,index_y));
        }
        else if (t1 == 0 && t2 == 1) {
            Index1D index_x(j1,k1,XBSpline);
            Index1D index_y(j2,k2,XWavelet);
            Lambda.insert(Index2D(index_x,index_y));
        }
        else if (t1 == 0 && t2 == 0) {
            Index1D index_x(j1,k1,XBSpline);
            Index1D index_y(j2,k2,XBSpline);
            Lambda.insert(Index2D(index_x,index_y));
        }
        else {
            std::cerr << "Could not read file." << std::endl;
            exit(1); return;
        }
    }
}

