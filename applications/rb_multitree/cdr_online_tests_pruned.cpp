#include "cdr_problem.h"
#include <fstream>
#include <sstream>

void
read_paramfile(string paramfilename, vector<ParamType>& v){
    T mu1, mu2;
    ifstream paramfile(paramfilename);
    if(paramfile.is_open()){
    	while(!paramfile.eof()){
        	paramfile >> mu1 >> mu2;
        	ParamType mu = {{mu1, mu2}};
        	v.push_back(mu);
    	}
    	paramfile.close();
    }
    else{
    	cerr << "Couldn't open " << paramfile << " for reading!" << endl;
    }
}

int main (int argc, char* argv[]) {

	//===============================================================//
	//========= PROBLEM SETUP  =======================//
	//===============================================================//
	
    int d   = 2;
    int d_  = 2;
    int j0  = 2;

    if(argc < 5 || argc > 6){
    	cerr << "Usage: " << argv[0] << " offline_data_folder snapshot_param_file output_name readSols(0/1) (Solution Folder)" << endl;
    	exit(1);
     }

    string offline_folder = argv[1];
    string param_file = argv[2];
    string output = argv[3];
    int readSols = atoi(argv[4]);
    string solfolder;
    if(readSols){
        solfolder = argv[5];
    }

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
    rb_system.read_alpha("S-5-5-Per-Intbc_Alpha.txt");

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
    rb_system.read_rb_data(offline_folder);
    rb_base.read_basisfunctions(offline_folder + "/bf");
    //rb_base.read_rieszrepresentors(offline_folder + "/representors");
    //rb_base.read_greedy_info("./Runs/Stage8/Run3/awgm_stage8_greedy_info.txt");

    // Read Test Parameters
    std::vector<ParamType> Xi_test, snapshot_params;
    
    read_paramfile("Xitest.txt", Xi_test);
    read_paramfile(param_file, snapshot_params);

    cout << "===============================" << endl;
    cout << "Test Parameters: " << endl;
    for(auto& mu : Xi_test){
    	cout << mu[0] << " " << mu[1] << endl;
    }
    cout << "===============================" << endl;

    // Test Reduced Solutions
    std::map<ParamType, std::vector<T> > errs, errbounds;

    
    // Construct pruned index sets
    int Nmax = snapshot_params.size();
    vector<vector<size_t> > pruned_indizes;
    
    DenseVectorT av_err(Nmax);
    DenseVectorT av_errest(Nmax);
    
    vector<size_t> pr_inds;
    pruned_indizes.push_back(pr_inds);
    for(size_t N = 1; N <= Nmax; N++){
        // Construct index set of pruned basis
        ParamType bf_param = snapshot_params[N-1];
        int index_first_occurence = -1;
        for(size_t i = 0; i < pruned_indizes[N-1].size(); ++i){
            if(fabs(snapshot_params[pruned_indizes[N-1][i]][0] - bf_param[0]) <= 1e-05 
            && fabs(snapshot_params[pruned_indizes[N-1][i]][1] - bf_param[1]) <= 1e-05){
                index_first_occurence = i;
                break;
            }
        }
        if(index_first_occurence >= 0){
            pruned_indizes[N-1][index_first_occurence] = N-1;
        }
        else{
            pruned_indizes[N-1].push_back(N-1);
        }
        
        cout << "Pruned basis has indices { ";
        for(auto& ind : pruned_indizes[N-1]){
            cout << ind << " ";
        }
        cout << "}" << endl << endl;
        
        pruned_indizes.push_back(pruned_indizes[N-1]);
    }

    DenseVectorT av_alt(3);
    // For all parameters:
    int count = 0;
    for(auto& mu : Xi_test){
        count++;
        
    	stringstream ufilename;
    	if(argc == 6){
    		ufilename << solfolder << "/";
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
        cout << " Mu = " << mu[0] << " " << mu[1] << endl;
        
    	for(size_t n = 1; n <= Nmax; ++n){

     		cout << n << endl << endl;
    		// Compute RB solution
    		DenseVectorT u_N = rb_system.get_rb_solution(pruned_indizes[n-1], mu);
            //cout << "u_" << n << " = " << u_N;
            
            /*// Test
            vector<size_t> rev_indices;
            for(int i = pruned_indizes[n-1].size()-1; i >= 0; --i){
                rev_indices.push_back(pruned_indizes[n-1][i]);
            }
            cout << "Reversed basis has indices { ";
            for(auto& ind : rev_indices){
                cout << ind << " ";
            }
            cout << "}" << endl << endl;
            DenseVectorT u_N_rev = rb_system.get_rb_solution(rev_indices, mu);
    		cout << "u_" << n << "_rev = " << u_N_rev; */
            

    		// Compute real error
    		DataType u_approx = rb_base.reconstruct_u_N(u_N, pruned_indizes[n-1]);
    		u_approx -= u;
    		err = u_approx.norm(2.);
    		errs_mu.push_back(err);
    		cout << "Err = " << u_approx.norm(2.)  << endl;
    		
    		if(n == 5){
        		DenseVectorT u_N_full = rb_system.get_rb_solution(n, mu);
                cout << "u_" << n << "_alt = " << u_N_full;
                
                DataType u_appr_full = rb_base.reconstruct_u_N(u_N_full,n);
                u_appr_full -= u;
                cout << "Err = " << u_appr_full.norm(2.) << endl;
                
                av_alt(1) += u_appr_full.norm(2.);
    		}
    		if(n == 6){
                vector<size_t> indizes = {{0,1,2,3,5}};
    		    DenseVectorT u_N_full = rb_system.get_rb_solution(indizes, mu);
                cout << "u_" << n << "_alt = " << u_N_full;
                
                DataType u_appr_full = rb_base.reconstruct_u_N(u_N_full,indizes);
                u_appr_full -= u;
                cout << "Err = " << u_appr_full.norm(2.) << endl;
                
                av_alt(2) += u_appr_full.norm(2.);
    		}
    		if(n == 7){
                vector<size_t> indizes = {{0,1,6,3,5}};
    		    DenseVectorT u_N_full = rb_system.get_rb_solution(indizes, mu);
                cout << "u_" << n << "_alt = " << u_N_full;
                
                DataType u_appr_full = rb_base.reconstruct_u_N(u_N_full,indizes);
                u_appr_full -= u;
                cout << "Err = " << u_appr_full.norm(2.) << endl;
                
                av_alt(3) += u_appr_full.norm(2.);
    		}
             
            av_err(n) += err;
            
    		// Compute error estimator
    		//errbound = rb_system.get_errorbound(pruned_indizes[n-1], u_N, mu);
    		errbounds_mu.push_back(errbound);
    		
            av_errest(n) += errbound;
    	}
    	errs.insert(make_pair(mu, errs_mu));
    	errbounds.insert(make_pair(mu, errbounds_mu));
    	    	
    }
    
    av_err *= 1./Xi_test.size();
    av_errest *= 1./Xi_test.size();
    av_alt *= 1./Xi_test.size();
    
   // rb_system.write_alpha("Test_Alphas.txt");

    ofstream errfile(output + "_PrunedOnlineTest_Errors.txt");
    if(errfile.is_open()){
        errfile << "# Mu1 Mu2";
        for(int n = 0; n < Nmax; ++n){
            errfile << " " << pruned_indizes[n].size();
        }
        errfile << endl;
    	
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

    ofstream errestfile(output + "_PrunedOnlineTest_ErrorEstimators.txt");
    if(errestfile.is_open()){
        errestfile << "# Mu1 Mu2";
        for(int n = 0; n < Nmax; ++n){
            errestfile << " " << pruned_indizes[n].size();
        }
        errestfile << endl;    	
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
    
    string av_file_name = output;
    av_file_name += "_PrunedOnlineTest_Averages.txt";
    ofstream av_file(av_file_name.c_str());
    if(av_file.is_open()){
        av_file << "# N_pruned avErr avErrEst av_err_alt" << endl;
        for(int i = 1; i <= Nmax; ++i){
            av_file << i << " " << pruned_indizes[i-1].size() << " " << av_err(i) << " " << av_errest(i);
            if(i>=5){
                av_file << " " << av_alt(i-4) << endl;
            }
            else{
                av_file << endl;
            }
        }
        
        av_file.close();        
    }
    else{
        cerr << "Couldn't open file " << av_file_name << " for writing! " << endl;
    }
}

