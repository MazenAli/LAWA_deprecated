#include "cdr_problem.h"
#include <fstream>
#include <sstream>

int main (int argc, char* argv[]) {

	//===============================================================//
	//========= PROBLEM SETUP  =======================//
	//===============================================================//
	
    int d   = 2;
    int d_  = 2;
    int j0  = 2;

    if(argc < 2 || argc > 3){
    	cerr << "Usage: " << argv[0] << " readSols(0/1) (Solution Folder)" << endl;
    	exit(1);
     }

    int readSols = atoi(argv[1]);

    //getchar();

    /// Basis initialization
    TrialBasis_Time      basis_per(d,d_,j0);
    TestBasis_Time       basis_int(d,d_,j0);
    Basis_Space 		 basis_intbc(d,d_,j0);
    basis_intbc.enforceBoundaryCondition<DirichletBC>();

    Basis2D_Trial basis2d_trial(basis_per,basis_intbc);
    Basis2D_Test  basis2d_test(basis_int,basis_intbc);

    /// Initialization of operators
    DenseVectorT no_singPts;
    Function<T> zero_Fct(zero_fct,no_singPts);
    Function<T> w_conv_Fct(weight_convection,no_singPts);

    // Bilinear Forms
    Convection1D_Time			ConvectionBil_t(basis_per, basis_int);
    Identity1D_Time 		    IdentityBil_t(basis_per, basis_int);
    Identity1D_Space 	        IdentityBil_x(basis_intbc, basis_intbc);
    Laplace1D_Space 	        LaplaceBil_x(basis_intbc, basis_intbc);
    Convection1D_Space 	        ConvectionBil_x(basis_intbc, basis_intbc,
    											zero_Fct, w_conv_Fct, zero_Fct, 10,true, false, true);
    Identity1D_Space 	        IdentityBil_Test_t(basis_int, basis_int);
    Identity1D_Time_Trial		IdentityBil_Trial_t(basis_per, basis_per);
    Convection1D_Time_Trial		ConvectionBil_Trial_t(basis_per, basis_per);


    RefConvection1D_Time 		RefConvectionBil_t(basis_per.refinementbasis, basis_int.refinementbasis);
    RefIdentity1D_Time 		    RefIdentityBil_t(basis_per.refinementbasis, basis_int.refinementbasis);
    RefIdentity1D_Space 	    RefIdentityBil_x(basis_intbc.refinementbasis, basis_intbc.refinementbasis);
    RefLaplace1D_Space 	        RefLaplaceBil_x(basis_intbc.refinementbasis, basis_intbc.refinementbasis);
    RefConvection1D_Space 	    RefConvectionBil_x(basis_intbc.refinementbasis, basis_intbc.refinementbasis,
													zero_Fct, w_conv_Fct, zero_Fct, 10,true, false, true);
    RefIdentity1D_Space 	    RefIdentityBil_Test_t(basis_int.refinementbasis, basis_int.refinementbasis);
    RefIdentity1D_Time_Trial	RefIdentityBil_Trial_t(basis_per.refinementbasis, basis_per.refinementbasis);
    RefConvection1D_Time_Trial	RefConvectionBil_Trial_t(basis_per.refinementbasis, basis_per.refinementbasis);

    // Transposed Bilinear Forms
    TranspConvection1D_Time 	TranspConvectionBil_t(basis_per, basis_int);
    TranspIdentity1D_Time 		TranspIdentityBil_t(basis_per, basis_int);
    TranspIdentity1D_Space 	    TranspIdentityBil_x(basis_intbc, basis_intbc);
    TranspLaplace1D_Space 	    TranspLaplaceBil_x(basis_intbc, basis_intbc);
    TranspConvection1D_Space 	TranspConvectionBil_x(basis_intbc, basis_intbc,
    													zero_Fct, w_conv_Fct, zero_Fct, 10,true, false, true);
    
    RefTranspConvection1D_Time 	RefTranspConvectionBil_t(basis_per.refinementbasis, basis_int.refinementbasis);
    RefTranspIdentity1D_Time 	RefTranspIdentityBil_t(basis_per.refinementbasis, basis_int.refinementbasis);
    RefTranspIdentity1D_Space 	RefTranspIdentityBil_x(basis_intbc.refinementbasis, basis_intbc.refinementbasis);
    RefTranspLaplace1D_Space 	RefTranspLaplaceBil_x(basis_intbc.refinementbasis, basis_intbc.refinementbasis);
    RefTranspConvection1D_Space RefTranspConvectionBil_x(basis_intbc.refinementbasis, basis_intbc.refinementbasis,
    													zero_Fct, w_conv_Fct, zero_Fct, 10, true, false, true);

    /// Initialization of local operator
    LOp_Conv1D_Time			lOp_Conv1D_t	(basis_int, basis_per, RefConvectionBil_t, ConvectionBil_t);
    LOp_Id1D_Time			lOp_Id1D_t  	(basis_int, basis_per, RefIdentityBil_t, IdentityBil_t);
    LOp_Id1D_Space			lOp_Id1D_x  	(basis_intbc, basis_intbc, RefIdentityBil_x, IdentityBil_x);
    LOp_Lapl1D_Space		lOp_Lapl1D_x	(basis_intbc, basis_intbc, RefLaplaceBil_x, LaplaceBil_x);
    LOp_Conv1D_Space		lOp_Conv1D_x	(basis_intbc, basis_intbc, RefConvectionBil_x, ConvectionBil_x);
    LOp_Id1D_Space			lOp_Id1D_Test_t (basis_int, basis_int, RefIdentityBil_Test_t, IdentityBil_Test_t);
    LOp_Id1D_Time_Trial		lOp_Id1D_Trial_t(basis_per, basis_per, RefIdentityBil_Trial_t, IdentityBil_Trial_t);
    LOp_Conv1D_Time_Trial 	lOp_Conv1D_Trial_t(basis_per, basis_per, RefConvectionBil_Trial_t, ConvectionBil_Trial_t);
    
    LOpT_Conv1D_Time		lOpT_Conv1D_t(basis_per, basis_int, RefTranspConvectionBil_t, TranspConvectionBil_t);
    LOpT_Id1D_Time			lOpT_Id1D_t  (basis_per, basis_int, RefTranspIdentityBil_t, TranspIdentityBil_t);
    LOpT_Id1D_Space			lOpT_Id1D_x  (basis_intbc, basis_intbc, RefTranspIdentityBil_x, TranspIdentityBil_x);
    LOpT_Lapl1D_Space		lOpT_Lapl1D_x(basis_intbc, basis_intbc, RefTranspLaplaceBil_x, TranspLaplaceBil_x);
    LOpT_Conv1D_Space		lOpT_Conv1D_x(basis_intbc, basis_intbc, RefTranspConvectionBil_x, TranspConvectionBil_x);

    LOp_Conv_Id_2D			localConvectionIdentityOp2D		(lOp_Conv1D_t, 		lOp_Id1D_x);
    LOp_Id_Id_2D			localIdentityIdentityOp2D		(lOp_Id1D_t, 		lOp_Id1D_x);
    LOp_Id_Lapl_2D			localIdentityLaplaceOp2D		(lOp_Id1D_t, 		lOp_Lapl1D_x);
    LOp_Id_Conv_2D			localIdentityConvectionOp2D		(lOp_Id1D_t, 		lOp_Conv1D_x);
    LOp_Test_Id_Id_2D		localIdentityIdentityOp2D_Test	(lOp_Id1D_Test_t, 	lOp_Id1D_x);
    LOp_Test_Id_Lapl_2D		localIdentityLaplaceOp2D_Test	(lOp_Id1D_Test_t, 	lOp_Lapl1D_x);
    LOp_Trial_Id_Id_2D		localIdentityIdentityOp2D_Trial	(lOp_Id1D_Trial_t, 	lOp_Id1D_x);
    LOp_Trial_Id_Lapl_2D 	localIdentityLaplaceOp2D_Trial	(lOp_Id1D_Trial_t, 	lOp_Lapl1D_x);
    LOp_Trial_Id_Conv_2D 	localIdentityConvectionOp2D_Trial(lOp_Id1D_Trial_t, lOp_Conv1D_x);
    LOp_Trial_Conv_Id_2D 	localConvectionIdentityOp2D_Trial(lOp_Conv1D_Trial_t,lOp_Id1D_x);

    LOpT_Conv_Id_2D		transpLocalConvectionIdentityOp2D	(lOpT_Conv1D_t, lOpT_Id1D_x);
    LOpT_Id_Id_2D		transpLocalIdentityIdentityOp2D		(lOpT_Id1D_t, 	lOpT_Id1D_x);
    LOpT_Id_Lapl_2D		transpLocalIdentityLaplaceOp2D		(lOpT_Id1D_t, 	lOpT_Lapl1D_x);
    LOpT_Id_Conv_2D		transpLocalIdentityConvectionOp2D	(lOpT_Id1D_t, 	lOpT_Conv1D_x);

    localConvectionIdentityOp2D.setJ(9);
    localIdentityIdentityOp2D.setJ(9);
    localIdentityLaplaceOp2D.setJ(9);
    localIdentityConvectionOp2D.setJ(9);
    transpLocalConvectionIdentityOp2D.setJ(9);
    transpLocalIdentityIdentityOp2D.setJ(9);
    transpLocalIdentityLaplaceOp2D.setJ(9);
    transpLocalIdentityConvectionOp2D.setJ(9);
    localIdentityIdentityOp2D_Test.setJ(9);
    localIdentityLaplaceOp2D_Test.setJ(9);
    localIdentityIdentityOp2D_Trial.setJ(9);
    localIdentityConvectionOp2D_Trial.setJ(9);
    localIdentityLaplaceOp2D_Trial.setJ(9);
    localConvectionIdentityOp2D_Trial.setJ(9);

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

    SeparableRhsIntegral2D_Trial	rhs_u(basis2d_trial, F1, nodeltas, nodeltas, 20);
    SeparableRhs_Trial           	F_u(rhs_u,noPrec);


	//===============================================================//
	//===============  RB SETUP =====================================//
	//===============================================================//

    // Affine Decompositions:
    // 	Left Hand Side
    vector<ThetaStructure<ParamType>::ThetaFct>	lhs_theta_fcts;
    lhs_theta_fcts.push_back(no_theta);
    lhs_theta_fcts.push_back(no_theta);
    lhs_theta_fcts.push_back(theta_conv);
    lhs_theta_fcts.push_back(theta_reac);
    ThetaStructure<ParamType> lhs_theta(lhs_theta_fcts);

    vector<AbstractLocalOperator2D<T>* > lhs_ops, lhs_opsT;
    lhs_ops.push_back(&localConvectionIdentityOp2D);
    lhs_ops.push_back(&localIdentityLaplaceOp2D);
    lhs_ops.push_back(&localIdentityConvectionOp2D);
    lhs_ops.push_back(&localIdentityIdentityOp2D);
    lhs_opsT.push_back(&transpLocalConvectionIdentityOp2D);
    lhs_opsT.push_back(&transpLocalIdentityLaplaceOp2D);
    lhs_opsT.push_back(&transpLocalIdentityConvectionOp2D);
    lhs_opsT.push_back(&transpLocalIdentityIdentityOp2D);

    Affine_Op_2D affine_lhs(lhs_theta, lhs_ops);
    Affine_Op_2D affine_lhs_T(lhs_theta, lhs_opsT);

    COp_TestInnProd innprod_Y	 (localIdentityIdentityOp2D_Test, localIdentityLaplaceOp2D_Test);
    TrialInnProdY 	innprod_Y_u_u(localIdentityIdentityOp2D_Trial, localIdentityLaplaceOp2D_Trial);

    vector<AbstractLocalOperator2D<T>* > A_u_u_ops;
    A_u_u_ops.push_back(&localConvectionIdentityOp2D_Trial);
    A_u_u_ops.push_back(&localIdentityLaplaceOp2D_Trial);
    A_u_u_ops.push_back(&localIdentityConvectionOp2D_Trial);
    A_u_u_ops.push_back(&localIdentityIdentityOp2D_Trial);
    Flex_COp_2D A_u_u(A_u_u_ops);

    // Right Hand Side
    vector<ThetaStructure<ParamType>::ThetaFct> rhs_theta_fcts;
    rhs_theta_fcts.push_back(no_theta);
    ThetaStructure<ParamType> rhs_theta(rhs_theta_fcts);
    vector<SeparableRhs*> rhs_fcts;
    rhs_fcts.push_back(&F);

    Affine_Rhs_2D affine_rhs(rhs_theta, rhs_fcts);
    RieszF_Rhs_2D rieszF_rhs(rhs_fcts);
    RieszA_Rhs_2D rieszA_rhs(lhs_ops);

    vector<SeparableRhs_Trial*> rhs_fcts_u;
    rhs_fcts_u.push_back(&F_u);
    Flex_Rhs_2D flex_rhs_u(rhs_fcts_u);


	//===============================================================//
	//===============  AWGM =========================================//
	//===============================================================//

    /* AWGM PG Parameters Default Values
    double tol = 5e-03;
	double alpha = 0.7;
	size_t max_its = 100;
	size_t max_basissize = 400000;
	bool reset_res = false;
	StableExpansionVersion stable_exp_u = FullExpansion,
    StableExpansionVersion stable_exp_res = FullExpansion,
    ResidualConstruction res_construction = SimpleStableExpansion,
	bool print_info = true;
	bool verbose = true;
	bool plot_solution = false;
	bool verbose_extra = false; //(print added wavelet indizes)
	size_t hashmapsize_trial = 10;
	size_t hashmapsize_test = 10;
	std::string info_filename = "awgm_cgls_conv_info.txt";
	std::string plot_filename = "awgm_cgls_u_plot";
	bool write_intermediary_solutions = false,
    std::string intermediary_solutions_filename = "awgm_cgls_u"
	*/

    /* IS Parameters Default Values
	bool adaptive_tol = true;
	size_t max_its = 100;
	double init_tol = 0.001;
	double res_reduction = 0.01;
	double absolute_tol = 1e-8;
	bool verbose = true;
	*/

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

    AWGM_PG_Parameters awgm_truth_parameters;
    IS_Parameters is_parameters;
    AWGM_Parameters awgm_riesz_f_parameters, awgm_riesz_a_parameters;
    awgm_truth_parameters.max_its = 1000;
    
    is_parameters.adaptive_tol = false;
    is_parameters.absolute_tol = 1e-08;


    //----------- Solver ---------------- //

    T gamma = 0.2;
    IndexSet<Index2D> LambdaTrial, LambdaTest;
    getSparseGridIndexSet(basis2d_trial,LambdaTrial,2,0,gamma);
    getSparseGridIndexSet(basis2d_test ,LambdaTest ,2,1,gamma);

    MT_AWGM_Truth awgm_u(basis2d_trial, basis2d_test, affine_lhs, affine_lhs_T,
   							 affine_rhs, rightPrec, leftPrec, awgm_truth_parameters, is_parameters);
    awgm_u.set_sol(dummy);
    awgm_u.awgm_params.tol = 1e-04;
    awgm_u.set_initial_indexsets(LambdaTrial,LambdaTest);


    MT_AWGM_Riesz_F awgm_rieszF(basis2d_test, innprod_Y, rieszF_rhs, leftPrec, awgm_riesz_f_parameters, is_parameters);
    awgm_rieszF.set_sol(dummy);
    awgm_rieszF.set_initial_indexset(LambdaTest);
    awgm_rieszF.awgm_params.tol = 5e-04;
    awgm_rieszF.awgm_params.info_filename = "awgm_stage9_rieszF_conv_info.txt";


    MT_AWGM_Riesz_A awgm_rieszA(basis2d_test, innprod_Y, rieszA_rhs, leftPrec, awgm_riesz_a_parameters, is_parameters);
    awgm_rieszA.set_sol(dummy);
    awgm_rieszA.set_initial_indexset(LambdaTest);
    awgm_rieszA.awgm_params.tol = 5e-04;
    awgm_rieszA.awgm_params.info_filename = "awgm_stage9_rieszA_conv_info.txt";

    MTTruthSolver rb_truth(awgm_u, awgm_rieszF, awgm_rieszA, innprod_Y_u_u, A_u_u, flex_rhs_u);


    //----------- RB System ---------------- //


    LB_Base<ParamType, MTTruthSolver> lb_base(rb_truth, lhs_theta);
    IndexSet<Index2D> LambdaTrial_Alpha_sparse, LambdaTrial_Alpha_full ;

    getSparseGridIndexSet(basis2d_trial,LambdaTrial_Alpha_sparse,3,0,gamma);
    //getFullIndexSet(basis2d_trial, LambdaTrial_Alpha_full, 3,3,0);


    cout << "++ Assembling Matrices for Alpha Computation ... " << endl << endl;
    lb_base.assemble_matrices_for_alpha_computation(LambdaTrial_Alpha_sparse);

    RB_Model rb_system(lhs_theta, rhs_theta, lb_base);
    //rb_system.read_alpha("S-5-5-Per-Intbc_Alpha.txt");

    rb_system.rb_params.ref_param = {{1., 1.}};
    rb_system.rb_params.call = call_gmres;

    //----------- RB Base ---------------- //

    RB_BaseModel rb_base(rb_system, rb_truth);

    /* RB Greedy Parameters Default Values
		TrainingType training_type = weak,
		double tol = 1e-2,
		std::size_t Nmax = 20,
		ParamType_min_param = ParamType(),
		ParamType max_param = ParamType(),
		intArray  training_params_per_dim = intArray(),
		bool print_info = true,
		std::string print_file = "greedy_info.txt",
		bool verbose = true,
		bool write_during_training = true,
		std::string trainingdata_folder = "training_data",
		bool print_paramset = false,
		bool erase_snapshot_params = false,
		bool orthonormalize_bfs = true,
		bool tighten_tol	= false,
		SnapshotTolReductionCrit snapshot_tol_red_crit = repeated_param,
		bool tighten_tol_rieszA = false,
		bool tighten_tol_rieszF = false,
		double tighten_tol_reduction = 0.1,
		bool update_snapshot = false,
		bool update_rieszF = false,
		bool update_rieszA = false,
		bool coarsen_rieszA_for_update = false,
		bool test_estimator_equivalence = false,
		bool equivalence_tol_factor = 1.,
		bool write_direct_representors = false,
		T min_error_reduction = 0.5);
     */


    /* RB Parameters Default Values
      		SolverCall call = call_cg,
			ParamType ref_param = ParamType(),
			bool verbose = true
     */

    ParamType mu_min = {{0., -9.}};
    ParamType mu_max = {{30, 15}};

    rb_base.greedy_params.tol = 1e-04;
    rb_base.greedy_params.min_param = mu_min;
    rb_base.greedy_params.max_param = mu_max;
    rb_base.greedy_params.Nmax = 	15;
    rb_base.greedy_params.nb_training_params = {{20, 20}};
    rb_base.greedy_params.print_paramset = true;
    rb_base.greedy_params.erase_snapshot_params = false;
    rb_base.greedy_params.orthonormalize_bfs = false;
    rb_base.greedy_params.print_file = "awgm_stage9_greedy_info.txt";
    rb_base.greedy_params.trainingdata_folder = "training_data_stage9";
    rb_base.greedy_params.tighten_tol = true;
    rb_base.greedy_params.snapshot_tol_red_crit = conv_rate_degradation;
    rb_base.greedy_params.tighten_tol_rieszA = false;
    rb_base.greedy_params.tighten_tol_rieszF = false;
    rb_base.greedy_params.update_snapshot = true;
    cout << "Parameters Training: " << std::endl << std::endl;
    rb_base.greedy_params.print();
    rb_system.rb_params.print();

    cout << "Parameters Truth Solver: " << std::endl << std::endl;
    awgm_u.awgm_params.print();

    awgm_u.is_params.print();

    cout << "Parameters Riesz Solver F : " << std::endl << std::endl;
    awgm_rieszF.awgm_params.print();

    cout << "Parameters Riesz Solver A : " << std::endl << std::endl;
    awgm_rieszA.awgm_params.print();
    
    // Read RB Data
    rb_system.read_rb_data("./Runs/Stage6/Run2/training_data_stage6");
    rb_base.read_basisfunctions("./Runs/Stage6/Run2/training_data_stage6/bf");
    rb_base.read_rieszrepresentors("./Runs/Stage6/Run2/training_data_stage6/representors");
    rb_base.read_greedy_info("./Runs/Stage6/Run2/awgm_stage6_greedy_info.txt");

    // Read Test Parameters
    std::vector<ParamType> Xi_test;

    ifstream paramfile("Xitest.txt");
    T mu1, mu2;
    if(paramfile.is_open()){
    	while(!paramfile.eof()){
        	paramfile >> mu1 >> mu2;
        	ParamType mu = {{mu1, mu2}};
        	Xi_test.push_back(mu);
    	}
    	paramfile.close();
    }
    else{
    	cerr << "Couldn't open Xitest.txt for reading!" << endl;
    }

    cout << "===============================" << endl;
    cout << "Test Parameters: " << endl;
    for(auto& mu : Xi_test){
    	cout << mu[0] << " " << mu[1] << endl;
    }
    cout << "===============================" << endl;


    // Test Reduced Solutions
    std::map<ParamType, std::vector<T> > errs, errbounds;

    // For all parameters:
    for(auto& mu : Xi_test){
        
    	stringstream ufilename;
    	if(argc == 3){
    		ufilename << argv[2] << "/";
    	}
    	ufilename << "u_" << mu[0] << "_" << mu[1] << ".txt";

    	DataType u;
    	if(!readSols){
        	// Compute truth solution
        	u = rb_truth.get_truth_solution(mu);
        	saveCoeffVector2D(u, basis2d_trial, ufilename.str().c_str());
    	}
    	else{
    		readCoeffVector2D(u, ufilename.str().c_str());
    	}


    	std::vector<T> errs_mu, errbounds_mu;

    	T err, errbound;
    	for(size_t n = 1; n <= (size_t)rb_system.RB_inner_product.numRows(); ++n){

    		// Compute RB solution
    		DenseVectorT u_N = rb_system.get_rb_solution(n, mu);

    		// Compute real error
    		DataType u_approx = rb_base.reconstruct_u_N(u_N, n);
    		u_approx -= u;
    		err = u_approx.norm(2.);
    		errs_mu.push_back(err);

    		// Compute error estimator
    		errbound = rb_system.get_errorbound(u_N, mu);
    		errbounds_mu.push_back(errbound);
    	}

    	errs.insert(make_pair(mu, errs_mu));
    	errbounds.insert(make_pair(mu, errbounds_mu));
    }
    
    rb_system.write_alpha("Test_Alphas.txt");

    ofstream errfile("Test_Errors.txt");
    if(errfile.is_open()){
    	errfile << "Mu1 Mu2 Err(N=1) Err(N=2) ...." << endl;
    	for(auto& entry : errs){
    		errfile << entry.first[0] << " " << entry.first[1];
    		for(auto& e : entry.second){
    			errfile << " " << e;
    		}
    		errfile << endl;
    	}
    	errfile.close();
    }
    else{
    	cerr << "Couldn't open Test_Errors.txt for writing!" << endl;
    }

    ofstream errestfile("Test_ErrorEstimators.txt");
    if(errestfile.is_open()){
    	errestfile << "Mu1 Mu2 ErrEst(N=1) ErrEst(N=2) ...." << endl;
    	for(auto& entry : errbounds){
    		errestfile << entry.first[0] << " " << entry.first[1];
    		for(auto& e : entry.second){
    			errestfile << " " << e;
    		}
    		errestfile << endl;
    	}
    	errestfile.close();
    }
    else{
    	cerr << "Couldn't open Test_Errors.txt for writing!" << endl;
    }
}

