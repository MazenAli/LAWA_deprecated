#include "cdr_problem.h"

template<typename _T, typename _IndexType>
class EmptyTruth{

public:

	typedef _T T;
	typedef _IndexType IndexType;

	_T lhs_u_u(std::size_t i, const IndexType& ind_row, const IndexType& ind_col) { return 0.;}

    _T innprod_Y_u_u(const IndexType& ind_row, const IndexType& ind_col) {return 0;}
};

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


	//===============================================================//
	//===============  AWGM =========================================//
	//===============================================================//



    //----------- RB System ---------------- //

    EmptyTruth<T,Index2D> empty_truth;

    LB_Base<ParamType, EmptyTruth<T,Index2D> > lb_base(empty_truth, lhs_theta);


    typedef Simple_RB_System<T,ParamType,  LB_Base<ParamType, EmptyTruth<T,Index2D> > > RBStandaloneSystem;
    RBStandaloneSystem rb_system(lhs_theta, rhs_theta, lb_base);

    rb_system.read_alpha("S-5-5-Per-Intbc_Alpha.txt");

    rb_system.rb_params.ref_param = {{1., 1.}};
    rb_system.rb_params.call = call_gmres;


    rb_system.read_rb_data("Runs/Stage6/Run2/training_data_stage6");

    ParamType mu = {{0,-9}};

    cout << "===== First N Basis Functions ====== " << endl << endl;

    for(int N = 4; N <= 6; N++){
    	cout << "-- N = " << N << " -- " << endl;
        DenseVectorT u_N_1 = rb_system.get_rb_solution(N, mu);
        cout << "   u_N = " << u_N_1 << endl;
        cout << "   ResDualNorm = " << rb_system.residual_dual_norm(u_N_1, mu) << endl << endl;
    }


    cout << "===== Selected Basis Functions ====== " << endl << endl;

    vector<size_t> indizes;
    indizes.push_back(0);
    indizes.push_back(1);
    indizes.push_back(2);
    indizes.push_back(3);

    //cout << "-- I = {1,2,3,4} -- " << endl;

    cout << "-- I = {";
    for(auto& i : indizes){
    	cout << i+1 << ",";
    }
    cout << "} -- " << endl;

    DenseVectorT u_N_2 = rb_system.get_rb_solution(indizes, mu);
    cout << "   u_N = " << u_N_2 << endl;
    cout << "   ResDualNorm = " << rb_system.residual_dual_norm(indizes, u_N_2, mu) << endl << endl;

    //cout << "-- I = {2,3,4,5} -- " << endl;

    indizes.erase(indizes.begin());
    indizes.push_back(4);
    cout << "-- I = {";
    for(auto& i : indizes){
    	cout << i+1 << ",";
    }
    cout << "} -- " << endl;


    u_N_2 = rb_system.get_rb_solution(indizes, mu);
    cout << "   u_N = " << u_N_2 << endl;
    cout << "   ResDualNorm = " << rb_system.residual_dual_norm(indizes, u_N_2, mu) << endl << endl;

    //cout << "-- I = {2,3,4,6} -- " << endl;

    indizes[3] = 5;
    cout << "-- I = {";
    for(auto& i : indizes){
    	cout << i+1 << ",";
    }
    cout << "} -- " << endl;
    u_N_2 = rb_system.get_rb_solution(indizes, mu);
    cout << "   u_N = " << u_N_2 << endl;
    cout << "   ResDualNorm = " << rb_system.residual_dual_norm(indizes, u_N_2, mu) << endl << endl;


    return 0;
}

