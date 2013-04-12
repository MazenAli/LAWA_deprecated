#include <iostream>
#include <iomanip>
#include <utility>
#include <array>
#include <vector>
#define _USE_MATH_DEFINES
#include <cmath>

#include <lawa/lawa.h>

/*
#include <lawa/methods/rb/datastructures/thetastructure.h>
#include <lawa/methods/rb/operators/affinelocaloperator.h>
#include <lawa/methods/rb/righthandsides/affinerhs.h>
#include <lawa/methods/adaptive/solvers/multitreeawgm2.h>
#include <lawa/methods/rb/righthandsides/flexiblebilformrhs.h>
#include <lawa/methods/rb/datastructures/mt_truth.h>
#include <lawa/methods/rb/datastructures/rb_system.h>
#include <lawa/methods/rb/datastructures/rb_base.h>
#include <lawa/methods/rb/datastructures/rb_parameters.h>
*/

#include <lawa/methods/rb/datastructures/lb_base.h>

using namespace std;
using namespace lawa;

//===============================================================//
//========= Simple RB System	  ===============================//
//===============================================================//

template <typename _T, typename _ParamType, typename LB_System>
class Simple_RB_System: public RB_System<_T,_ParamType> {

    typedef std::map<_ParamType, _T>   StabConstVec;

public:
	Simple_RB_System(ThetaStructure<_ParamType>& _thetas_a,
			  	  ThetaStructure<_ParamType>& _thetas_f,
			  	  LB_System& _lb_system)
	 : RB_System<_T,_ParamType>(_thetas_a, _thetas_f),
	   lb_system(_lb_system){}

    _T
	alpha_LB(_ParamType& mu)
    {
    	// If it is precalculated, return value
    	for(auto& el: alphas){
    	    bool is_mu = true;
			for(std::size_t i = 0; i < mu.size(); ++i){
				if( (mu[i]-(el.first)[i]) > 1e-4 ){
					is_mu = false;
					break;
				}
			}
			if(is_mu){
				return el.second;
			}
    	}
    	// Else we have to calculate it
    	_T alpha = lb_system.calculate_alpha(mu);
    	alphas.insert(std::make_pair(mu, alpha));
    	return alpha;
    }

    void
    read_alpha(const char* filename){
        ifstream alphaFile(filename);
        while(alphaFile.good()){
             _ParamType mu;
             _T alpha;
             size_t pdim = ParamInfo<_ParamType>::dim;
             for(size_t p = 0; p < pdim; ++p){
            	 alphaFile >> mu[p];
             }
             alphaFile >> alpha;
             alphas.insert(std::make_pair(mu, alpha));
        }
    }

    void
    write_alpha(const char* filename){
        size_t pdim = ParamInfo<_ParamType>::dim;
        ofstream alphaFile(filename);
        if(alphaFile.is_open()){
        	for(auto& a : alphas){
                for(size_t p = 0; p < pdim; ++p){
               	 alphaFile << a.first[p] << " ";
                }
                alphaFile << a.second << std::endl;
        	}
        }
        else{
        	std::cerr << "Couldn't open file " << filename << " for writing." << std::endl;
        	return;
        }
    }

private:
    StabConstVec	alphas;
    LB_System&		lb_system;
};


//===============================================================//
//========= TYPEDEFS  =======================//
//===============================================================//

//==== General ====//
typedef double T;
typedef flens::GeMatrix<flens::FullStorage<T, cxxblas::ColMajor> >  FullColMatrixT;
typedef flens::DenseVector<flens::Array<T> >                        DenseVectorT;
typedef Coefficients<Lexicographical,T,Index2D>						CoeffVector;

//==== Basis 1D & 2D ====//
typedef Basis<T, Primal, Periodic, CDF>	                            TrialBasis_Time;
//typedef Basis<T, Primal, Periodic, CDF>								TestBasis_Time;
typedef Basis<T, Primal, Interval, Dijkema>						    TestBasis_Time;
typedef Basis<T, Primal, Interval, Dijkema>							Basis_Space;
typedef TrialBasis_Time::RefinementBasis                            TrialRefBasis_Time;
typedef TestBasis_Time::RefinementBasis                           	TestRefBasis_Time;
typedef Basis_Space::RefinementBasis                                RefBasis_Space;

// !!! ATTENTION: Use typedefs in space for definition of space-time inner product !!!
//    => Assumption that test basis in time = space basis (+/- boundary conditions)

typedef TensorBasis2D<Adaptive,TrialBasis_Time,Basis_Space>         Basis2D_Trial;
typedef TensorBasis2D<Adaptive,TestBasis_Time,Basis_Space>          Basis2D_Test;

//==== Adaptive Operators ====//
typedef AdaptiveConvectionOperator1D_PG<T,TrialBasis_Time,TestBasis_Time>                   Convection1D_Time;
typedef AdaptiveIdentityOperator1D_PG<T,TrialBasis_Time,TestBasis_Time>                     Identity1D_Time;
typedef AdaptiveIdentityOperator1D_PG<T,Basis_Space,Basis_Space>                            Identity1D_Space;
typedef AdaptiveLaplaceOperator1D_PG<T,Basis_Space,Basis_Space>                             Laplace1D_Space;
typedef AdaptiveWeightedPDEOperator1D_PG<T,Basis_Space,Basis_Space>                         Convection1D_Space;
typedef AdaptiveIdentityOperator1D_PG<T,TrialBasis_Time,TrialBasis_Time>                    Identity1D_Time_Trial;
typedef AdaptiveConvectionOperator1D_PG<T,TrialBasis_Time,TrialBasis_Time>                  Convection1D_Time_Trial;

typedef TransposedAdaptiveConvectionOperator1D_PG<T,TrialBasis_Time,TestBasis_Time>         TranspConvection1D_Time;
typedef TransposedAdaptiveIdentityOperator1D_PG<T,TrialBasis_Time,TestBasis_Time>           TranspIdentity1D_Time;
typedef TransposedAdaptiveIdentityOperator1D_PG<T,Basis_Space,Basis_Space>                  TranspIdentity1D_Space;
typedef TransposedAdaptiveLaplaceOperator1D_PG<T,Basis_Space,Basis_Space>                   TranspLaplace1D_Space;
typedef TransposedAdaptiveWeightedPDEOperator1D_PG<T,Basis_Space,Basis_Space>               TranspConvection1D_Space;

typedef AdaptiveConvectionOperator1D_PG<T,TrialRefBasis_Time,TestRefBasis_Time>             RefConvection1D_Time;
typedef AdaptiveIdentityOperator1D_PG<T,TrialRefBasis_Time,TestRefBasis_Time>               RefIdentity1D_Time;
typedef AdaptiveIdentityOperator1D_PG<T,RefBasis_Space,RefBasis_Space>                      RefIdentity1D_Space;
typedef AdaptiveLaplaceOperator1D_PG<T,RefBasis_Space,RefBasis_Space>                       RefLaplace1D_Space;
typedef AdaptiveWeightedPDEOperator1D_PG<T,RefBasis_Space,RefBasis_Space>                   RefConvection1D_Space;
typedef AdaptiveIdentityOperator1D_PG<T,TrialRefBasis_Time,TrialRefBasis_Time>              RefIdentity1D_Time_Trial;
typedef AdaptiveConvectionOperator1D_PG<T,TrialRefBasis_Time,TrialRefBasis_Time>            RefConvection1D_Time_Trial;

typedef TransposedAdaptiveConvectionOperator1D_PG<T,TrialRefBasis_Time,TestRefBasis_Time>   RefTranspConvection1D_Time;
typedef TransposedAdaptiveIdentityOperator1D_PG<T,TrialRefBasis_Time,TestRefBasis_Time>     RefTranspIdentity1D_Time;
typedef TransposedAdaptiveIdentityOperator1D_PG<T,RefBasis_Space,RefBasis_Space>            RefTranspIdentity1D_Space;
typedef TransposedAdaptiveLaplaceOperator1D_PG<T,RefBasis_Space,RefBasis_Space>             RefTranspLaplace1D_Space;
typedef TransposedAdaptiveWeightedPDEOperator1D_PG<T,RefBasis_Space,RefBasis_Space>         RefTranspConvection1D_Space;

//==== LocalOperators ====//
typedef LocalOperator1D<TestBasis_Time,TrialBasis_Time, 
                        RefConvection1D_Time,Convection1D_Time>                 LOp_Conv1D_Time;
typedef LocalOperator1D<TestBasis_Time,TrialBasis_Time, 
                        RefIdentity1D_Time,Identity1D_Time>                     LOp_Id1D_Time;                        
typedef LocalOperator1D<Basis_Space,Basis_Space,
						RefIdentity1D_Space,Identity1D_Space>	                LOp_Id1D_Space;
typedef LocalOperator1D<Basis_Space,Basis_Space,
						RefLaplace1D_Space,Laplace1D_Space>	                    LOp_Lapl1D_Space;
typedef LocalOperator1D<Basis_Space,Basis_Space,
						RefConvection1D_Space,Convection1D_Space>	            LOp_Conv1D_Space;
typedef LocalOperator1D<TrialBasis_Time,TrialBasis_Time,
                        RefIdentity1D_Time_Trial,Identity1D_Time_Trial>         LOp_Id1D_Time_Trial;
typedef LocalOperator1D<TrialBasis_Time,TrialBasis_Time,
                        RefConvection1D_Time_Trial,Convection1D_Time_Trial>     LOp_Conv1D_Time_Trial;

typedef LocalOperator1D<TrialBasis_Time,TestBasis_Time, 
                        RefTranspConvection1D_Time,TranspConvection1D_Time>     LOpT_Conv1D_Time;
typedef LocalOperator1D<TrialBasis_Time,TestBasis_Time, 
                        RefTranspIdentity1D_Time,TranspIdentity1D_Time>         LOpT_Id1D_Time;                        
typedef LocalOperator1D<Basis_Space,Basis_Space,
						RefTranspIdentity1D_Space,TranspIdentity1D_Space>	    LOpT_Id1D_Space;
typedef LocalOperator1D<Basis_Space,Basis_Space,
						RefTranspLaplace1D_Space,TranspLaplace1D_Space>	        LOpT_Lapl1D_Space;
typedef LocalOperator1D<Basis_Space,Basis_Space,
						RefTranspConvection1D_Space,TranspConvection1D_Space>	LOpT_Conv1D_Space;

typedef LocalOperator2D<LOp_Conv1D_Time, LOp_Id1D_Space>                        LOp_Conv_Id_2D;
typedef LocalOperator2D<LOp_Id1D_Time, LOp_Id1D_Space>                        	LOp_Id_Id_2D;
typedef LocalOperator2D<LOp_Id1D_Time, LOp_Lapl1D_Space>                        LOp_Id_Lapl_2D;
typedef LocalOperator2D<LOp_Id1D_Time, LOp_Conv1D_Space>                        LOp_Id_Conv_2D;
typedef LocalOperator2D<LOpT_Conv1D_Time, LOpT_Id1D_Space>                      LOpT_Conv_Id_2D;
typedef LocalOperator2D<LOpT_Id1D_Time, LOpT_Id1D_Space>                      	LOpT_Id_Id_2D;
typedef LocalOperator2D<LOpT_Id1D_Time, LOpT_Lapl1D_Space>                      LOpT_Id_Lapl_2D;
typedef LocalOperator2D<LOpT_Id1D_Time, LOpT_Conv1D_Space>                      LOpT_Id_Conv_2D;


typedef LocalOperator2D<LOp_Id1D_Space,LOp_Id1D_Space>							LOp_Test_Id_Id_2D;
typedef LocalOperator2D<LOp_Id1D_Space,LOp_Lapl1D_Space>						LOp_Test_Id_Lapl_2D;
typedef LocalOperator2D<LOp_Id1D_Time_Trial,LOp_Id1D_Space>						LOp_Trial_Id_Id_2D;
typedef LocalOperator2D<LOp_Id1D_Time_Trial,LOp_Conv1D_Space>					LOp_Trial_Id_Conv_2D;
typedef LocalOperator2D<LOp_Id1D_Time_Trial,LOp_Lapl1D_Space>					LOp_Trial_Id_Lapl_2D;
typedef LocalOperator2D<LOp_Conv1D_Time_Trial,LOp_Id1D_Space>					LOp_Trial_Conv_Id_2D;

//==== CompoundOperators ====//
typedef FlexibleCompoundLocalOperator<Index2D,AbstractLocalOperator2D<T> > 		Flex_COp_2D;
typedef CompoundLocalOperator<Index2D,LOp_Conv_Id_2D,LOp_Id_Lapl_2D>    		COp_Heat;
typedef CompoundLocalOperator<Index2D,LOpT_Conv_Id_2D,LOpT_Id_Lapl_2D>    	    COpT_Heat;
typedef CompoundLocalOperator<Index2D,LOp_Test_Id_Id_2D, LOp_Test_Id_Lapl_2D>	COp_TestInnProd;
typedef CompoundLocalOperator<Index2D,LOp_Trial_Id_Id_2D, LOp_Trial_Id_Lapl_2D>	TrialInnProdY;

//==== Preconditioners ====//
typedef AdaptiveLeftNormPreconditioner2D<T,Basis2D_Test>            LeftPrec2D;
typedef AdaptiveRightNormPreconditioner2D_c<T,Basis2D_Trial>        RightPrec2D;
typedef NoPreconditioner<T, Index2D>								NoPrec2D;

//==== RightHandSides ====//
typedef SeparableRHS2D<T,Basis2D_Test>                              SeparableRhsIntegral2D;
typedef RHS<T,Index2D,SeparableRhsIntegral2D,
            NoPrec2D>                                         		SeparableRhs;
typedef SeparableRHS2D<T,Basis2D_Trial>                             SeparableRhsIntegral2D_Trial;
typedef RHS<T,Index2D,SeparableRhsIntegral2D_Trial,
            NoPrec2D>                                         		SeparableRhs_Trial;


//==== RB Stuff ====//
const size_t PDim = 2;

typedef CoeffVector 								 						DataType;
typedef array<T,PDim>	 													ParamType;

typedef AffineLocalOperator<Index2D,AbstractLocalOperator2D<T>,ParamType>	Affine_Op_2D;
typedef AffineRhs<T,Index2D,SeparableRhs,ParamType>							Affine_Rhs_2D;
typedef FlexibleCompoundRhs<T,Index2D,SeparableRhs>							RieszF_Rhs_2D;
typedef FlexibleBilformRhs<Index2D,AbstractLocalOperator2D<T> >				RieszA_Rhs_2D;
typedef FlexibleCompoundRhs<T,Index2D,SeparableRhs_Trial>					Flex_Rhs_2D;

typedef MultiTreeAWGM_PG<Index2D,Basis2D_Trial, Basis2D_Test,Affine_Op_2D,
		Affine_Op_2D,Affine_Rhs_2D,RightPrec2D,LeftPrec2D>					MT_AWGM_Truth;
typedef MultiTreeAWGM2<Index2D,Basis2D_Test, COp_TestInnProd,
		RieszF_Rhs_2D,LeftPrec2D>											MT_AWGM_Riesz_F;
typedef MultiTreeAWGM2<Index2D,Basis2D_Test, COp_TestInnProd,
		RieszA_Rhs_2D,LeftPrec2D>											MT_AWGM_Riesz_A;

typedef MT_Truth<DataType,ParamType,MT_AWGM_Truth,
				 MT_AWGM_Riesz_F,MT_AWGM_Riesz_A,TrialInnProdY,
				 Flex_COp_2D, Flex_Rhs_2D>									MTTruthSolver;

typedef LB_Base<ParamType, MTTruthSolver> 									LB_Type;
typedef Simple_RB_System<T,ParamType, LB_Type>								RB_Model;
typedef RB_Base<RB_Model,MTTruthSolver,DataType,ParamType>					RB_BaseModel;

T f_t(T t)
{
    return cos(2.L*M_PI*t);
}

/*T f_x(T x)
{
    return -4*(x-0.5)*(x-0.5) + 1;
}
*/
T f_x(T)
{
	return 1.;
}

T dummy(T, T)
{
    return 0;
}

T no_theta(const std::array<T,PDim>& /*mu*/)
{
	return 1.;
}

T theta_conv(const std::array<T,PDim>& mu)
{
	return mu[0];
}

T theta_reac(const std::array<T,PDim>& mu)
{
	return mu[1];
}

T
weight_convection(T x)
{
  return 0.5 - x;
}

T zero_fct(T /*x*/){
	return 0;
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
	bool print_info = true;
	bool verbose = true;
	bool plot_solution = false;
	bool verbose_extra = false; //(print added wavelet indizes)
	size_t hashmapsize_trial = 10;
	size_t hashmapsize_test = 10;
	std::string info_filename = "awgm_cgls_conv_info.txt";
	std::string plot_filename = "awgm_cgls_u_plot";
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
    awgm_u.awgm_params.tol = 1e-03;
    awgm_u.set_initial_indexsets(LambdaTrial,LambdaTest);


    MT_AWGM_Riesz_F awgm_rieszF(basis2d_test, innprod_Y, rieszF_rhs, leftPrec, awgm_riesz_f_parameters, is_parameters);
    awgm_rieszF.set_sol(dummy);
    awgm_rieszF.set_initial_indexset(LambdaTest);
    awgm_rieszF.awgm_params.tol = 5e-03;
    awgm_rieszF.awgm_params.info_filename = "awgm_stage4_rieszF_conv_info.txt";


    MT_AWGM_Riesz_A awgm_rieszA(basis2d_test, innprod_Y, rieszA_rhs, leftPrec, awgm_riesz_a_parameters, is_parameters);
    awgm_rieszA.set_sol(dummy);
    awgm_rieszA.set_initial_indexset(LambdaTest);
    awgm_rieszA.awgm_params.tol = 5e-03;
    awgm_rieszA.awgm_params.info_filename = "awgm_stage4_rieszA_conv_info.txt";

    MTTruthSolver rb_truth(awgm_u, awgm_rieszF, awgm_rieszA, innprod_Y_u_u, A_u_u, flex_rhs_u);


    //----------- RB System ---------------- //


    LB_Base<ParamType, MTTruthSolver> lb_base(rb_truth, lhs_theta);
    IndexSet<Index2D> LambdaTrial_Alpha_sparse, LambdaTrial_Alpha_full ;

    //getSparseGridIndexSet(basis2d_trial,LambdaTrial_Alpha_sparse,3,0,gamma);
    getFullIndexSet(basis2d_trial, LambdaTrial_Alpha_full, 3,3,0);


    cout << "++ Assembling Matrices for Alpha Computation ... " << endl << endl;
    lb_base.assemble_matrices_for_alpha_computation(LambdaTrial_Alpha_full);


    RB_Model rb_system(lhs_theta, rhs_theta, lb_base);

    rb_system.rb_params.ref_param = {{1., 1.}};
    rb_system.rb_params.call = call_gmres;

    //----------- RB Base ---------------- //

    RB_BaseModel rb_base(rb_system, rb_truth);

    /* RB Greedy Parameters Default Values
      		double tol = 1e-2,
			size_t Nmax = 20,
			ParamType min_param = ParamType(),
			ParamType max_param = ParamType(),
			intArray  training_params_per_dim = intArray(),
			bool print_info = true,
    		std::string print_file = "greedy_info.txt",
			bool verbose = true,
			bool write_during_training = true,
    		std::string trainingdata_folder = "training_data",
			bool print_paramset = false,
			bool erase_snapshot_params = false,
			bool orthonormalize_bfs = true
     */


    /* RB Parameters Default Values
      		SolverCall call = call_cg,
			ParamType ref_param = ParamType(),
			bool verbose = true
     */

    ParamType mu_min = {{0., -9.}};
    ParamType mu_max = {{30, 15}};

    rb_base.greedy_params.min_param = mu_min;
    rb_base.greedy_params.max_param = mu_max;
    rb_base.greedy_params.Nmax = 	1;
    rb_base.greedy_params.nb_training_params = {{20, 20}};
    rb_base.greedy_params.print_paramset = true;
    rb_base.greedy_params.erase_snapshot_params = false;
    rb_base.greedy_params.orthonormalize_bfs = false;
    rb_base.greedy_params.print_file = "awgm_greedy_info.txt";
    rb_base.greedy_params.trainingdata_folder = "training_data";

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

    rb_base.train_Greedy();

    rb_system.write_alpha("offline_data/alphas.txt");
    rb_system.write_rb_data("offline_data");
    rb_base.write_basisfunctions("offline_data");
    rb_base.write_rieszrepresentors("offline_data");


    return 0;
}

