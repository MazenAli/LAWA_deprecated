#include "cdr_problem.h"

template<typename _T, typename _IndexType>
class EmptyTruth{

public:

	typedef _T T;
	typedef _IndexType IndexType;

	_T lhs_u_u(std::size_t i, const IndexType& ind_row, const IndexType& ind_col) { return 0.;}

    _T innprod_Y_u_u(const IndexType& ind_row, const IndexType& ind_col) {return 0;}
};

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
	
	if(argc != 4){
        cerr << "Usage: " << argv[0] << " offline_data_folder snapshot_param_file output_name" << endl;
        exit(1);
	}
	
    string offline_folder = argv[1];
    string param_file = argv[2];
    string output = argv[3];
	
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

    /// Initialization of preconditioner
    LeftPrec2D leftPrec(basis2d_test);
    RightPrec2D rightPrec(basis2d_trial);

    NoPrec2D noPrec;


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

    // Right Hand Side
    vector<ThetaStructure<ParamType>::ThetaFct> rhs_theta_fcts;
    rhs_theta_fcts.push_back(no_theta);
    ThetaStructure<ParamType> rhs_theta(rhs_theta_fcts);

    //----------- RB System ---------------- //

    EmptyTruth<T,Index2D> empty_truth;

    LB_Base<ParamType, EmptyTruth<T,Index2D> > lb_base(empty_truth, lhs_theta);

    typedef Simple_RB_System<T,ParamType,  LB_Base<ParamType, EmptyTruth<T,Index2D> > > RBStandaloneSystem;
    RBStandaloneSystem rb_system(lhs_theta, rhs_theta, lb_base);

    rb_system.read_alpha("S-5-5-Per-Intbc_Alpha.txt");

    rb_system.rb_params.ref_param = {{1., 1.}};
    rb_system.rb_params.call = call_gmres;

    rb_system.read_rb_data(offline_folder);
    
    // Read Training Parameters
    std::vector<ParamType> Xi_train, snapshot_params;
    
    read_paramfile("Xitrain.txt", Xi_train);
    read_paramfile(param_file, snapshot_params);
    
    /*
    ifstream paramfile("Xitrain.txt");
    T mu1, mu2;
    if(paramfile.is_open()){
    	while(!paramfile.eof()){
        	paramfile >> mu1 >> mu2;
        	ParamType mu = {{mu1, mu2}};
        	Xi_train.push_back(mu);
    	}
    	paramfile.close();
    }
    else{
    	cerr << "Couldn't open Xitrain.txt for reading!" << endl;
    }
    
    // Read Snapshot Parameters
    ifstream paramfile("Xitrain.txt");
    T mu1, mu2;
    if(paramfile.is_open()){
    	while(!paramfile.eof()){
        	paramfile >> mu1 >> mu2;
        	ParamType mu = {{mu1, mu2}};
        	Xi_train.push_back(mu);
    	}
    	paramfile.close();
    }
    else{
    	cerr << "Couldn't open Xitrain.txt for reading!" << endl;
    }*/

	//===============================================================//
	//===============  ONLINE TESTS =================================//
	//===============================================================//
	
    int Nmax = rb_system.RB_inner_product.numRows();
    FullColMatrixT err_ests_full_basis(Xi_train.size(), Nmax);
    FullColMatrixT err_ests_pruned_basis(Xi_train.size() + 1, Nmax);
   
    DenseVectorT max_full_basis(Nmax);
    DenseVectorT max_pruned_basis(Nmax);

    vector<size_t> pruned_indizes;

    cout << "===== Testing the full basis ====== " << endl << endl;

    bool pruning = false;
    for(int N = 1; N <= Nmax; N++){
        cout << "--------------------------------------------" << endl;
    	cout << "--------- N = " << N << " ----------------- " << endl;
    	cout << "--------------------------------------------" << endl << endl;
    	
    	// Construct index set of pruned basis
        ParamType bf_param = snapshot_params[N-1];
        int index_first_occurence = -1;
        for(size_t i = 0; i < pruned_indizes.size(); ++i){
            if(fabs(snapshot_params[pruned_indizes[i]][0] - bf_param[0]) <= 1e-05 
            && fabs(snapshot_params[pruned_indizes[i]][1] - bf_param[1]) <= 1e-05){
                index_first_occurence = i;
                break;
            }
        }
        if(index_first_occurence >= 0){
            pruning = true;
            pruned_indizes[index_first_occurence] = N-1;
        }
        else{
            pruned_indizes.push_back(N-1);
        }
        cout << "Pruned basis has indices { ";
        for(auto& ind : pruned_indizes){
            cout << ind << " ";
        }
        cout << "}" << endl << endl;
        err_ests_pruned_basis(1, N) = pruned_indizes.size();
        
        int row_index = 1;
        for(auto& mu : Xi_train){
            cout << "  Mu = [" << mu[0] << ", " << mu[1] << "]" << endl;
            DenseVectorT u_N = rb_system.get_rb_solution(N, mu);  
            cout << "  u_N = " << u_N;
            err_ests_full_basis(row_index, N) = rb_system.get_errorbound(u_N, mu); 
            cout << "  Error Estimator: " << err_ests_full_basis(row_index, N) << endl;
            if(err_ests_full_basis(row_index,N) > max_full_basis(N)){
                max_full_basis(N) = err_ests_full_basis(row_index,N);
            }
            
            if(pruning){
                DenseVectorT u_N_pruned = rb_system.get_rb_solution(pruned_indizes, mu);
                cout << "  u_N_pruned = " << u_N_pruned;
                err_ests_pruned_basis(row_index+1, N) = rb_system.get_errorbound(pruned_indizes, u_N_pruned, mu); 
                cout << "  Error Estimator pruned: " << err_ests_pruned_basis(row_index+1, N) << endl;                
            }
            else{
                err_ests_pruned_basis(row_index+1, N) = err_ests_full_basis(row_index, N);
            }
            if(err_ests_pruned_basis(row_index+1, N) > max_pruned_basis(N)){
                max_pruned_basis(N) = err_ests_pruned_basis(row_index+1, N);
            }

            row_index++;         
        }
    }
      
    string full_file_name = output;
    full_file_name += "_fullBasis.txt";
    ofstream full_file(full_file_name.c_str());
    if(full_file.is_open()){
        full_file << err_ests_full_basis << endl;
        full_file.close();
    }
    else{
        cerr << "Couldn't open file " << full_file_name << " for writing! " << endl;
    }
    
    string pruned_file_name = output;
    pruned_file_name += "_prunedBasis.txt";
    ofstream pruned_file(pruned_file_name.c_str());
    if(pruned_file.is_open()){
        pruned_file << err_ests_pruned_basis << endl;
        pruned_file.close();
    }
    else{
        cerr << "Couldn't open file " << pruned_file_name << " for writing! " << endl;
    }

    string max_file_name = output;
    max_file_name += "_maxErrors.txt";
    ofstream max_file(max_file_name.c_str());
    if(max_file.is_open()){
        max_file << "# N_full max_full N_pruned max_pruned" << endl;
        for(int i = 1; i <= Nmax; ++i){
            max_file << i << " " << max_full_basis(i) << " " << err_ests_pruned_basis(1,i) << " " << max_pruned_basis(i) << endl;
        }
        
        max_file.close();        
    }
    else{
        cerr << "Couldn't open file " << max_file_name << " for writing! " << endl;
    }
    
    return 0;
}

