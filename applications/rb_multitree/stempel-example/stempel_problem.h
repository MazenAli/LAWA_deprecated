#include <iostream>
#include <iomanip>
#include <utility>
#include <array>
#include <vector>
#define _USE_MATH_DEFINES
#include <cmath>

#include <lawa/lawa.h>

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
    	/*// If it is precalculated, return value
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
    	*/
        return std::min(mu[0], 1.);
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
    
    void
    set_alpha(const StabConstVec& _alphas){
        alphas = _alphas;
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
const DomainType      DomainType_XY                                 = Interval;
const FunctionSide    FctSide_X                                     = Primal;
const Construction    Constr_X                                      = Dijkema;
const Construction    RefConstr_X                                   = Dijkema;
const FunctionSide    FctSide_Y                                     = Primal;
const Construction    Constr_Y                                      = Dijkema;
const Construction    RefConstr_Y                                   = Dijkema;
typedef Basis<T,FctSide_X, DomainType_XY, Constr_X>					Basis_X;
typedef Basis<T,FctSide_Y, DomainType_XY, Constr_Y>					Basis_Y;
typedef Basis_X::RefinementBasis                                    RefBasis_X;
typedef Basis_Y::RefinementBasis                           	        RefBasis_Y;

typedef TensorBasis2D<Adaptive,Basis_X,Basis_Y>                     Basis2D;

//==== Adaptive Operators ====//
typedef AdaptiveIdentityOperator1D<T,FctSide_X, DomainType_XY, Constr_X>                Identity1D_X;
typedef AdaptiveIdentityOperator1D<T,FctSide_Y, DomainType_XY, Constr_Y>                Identity1D_Y;
typedef AdaptiveLaplaceOperator1D<T,FctSide_X, DomainType_XY, Constr_X>                 Laplace1D_X;
typedef AdaptiveLaplaceOperator1D<T,FctSide_Y, DomainType_XY, Constr_Y>                 Laplace1D_Y;
typedef AdaptiveWeightedPDEOperator1D<T,FctSide_X, DomainType_XY, Constr_X>             WeightedLaplace1D_X;
typedef AdaptiveWeightedPDEOperator1D<T,FctSide_X, DomainType_XY, Constr_X>             WeightedIdentity1D_X;

typedef AdaptiveIdentityOperator1D<T,FctSide_X, DomainType_XY, RefConstr_X>             RefIdentity1D_X;
typedef AdaptiveIdentityOperator1D<T,FctSide_Y, DomainType_XY, RefConstr_Y>             RefIdentity1D_Y;
typedef AdaptiveLaplaceOperator1D<T,FctSide_X, DomainType_XY, RefConstr_X>              RefLaplace1D_X;
typedef AdaptiveLaplaceOperator1D<T,FctSide_Y, DomainType_XY, RefConstr_Y>              RefLaplace1D_Y;
typedef AdaptiveWeightedPDEOperator1D<T,FctSide_X, DomainType_XY, RefConstr_X>          RefWeightedLaplace1D_X;
typedef AdaptiveWeightedPDEOperator1D<T,FctSide_X, DomainType_XY, RefConstr_X>          RefWeightedIdentity1D_X;

//==== LocalOperators ====//
typedef LocalOperator1D<Basis_X,Basis_X,RefIdentity1D_X,Identity1D_X>                   LOp_Id1D_X;
typedef LocalOperator1D<Basis_Y,Basis_Y,RefIdentity1D_Y,Identity1D_Y>                   LOp_Id1D_Y;
typedef LocalOperator1D<Basis_X,Basis_X,RefLaplace1D_X,Laplace1D_X>                     LOp_Lapl1D_X;
typedef LocalOperator1D<Basis_Y,Basis_Y,RefLaplace1D_Y,Laplace1D_Y>                     LOp_Lapl1D_Y;
typedef LocalOperator1D<Basis_X,Basis_X,RefWeightedLaplace1D_X,WeightedLaplace1D_X>     LOp_WLapl1D_X;
typedef LocalOperator1D<Basis_X,Basis_X,RefWeightedIdentity1D_X,WeightedIdentity1D_X>   LOp_WId1D_X;

typedef LocalOperator2D<LOp_Id1D_X, LOp_Id1D_Y>                        	                LOp_Id_Id_2D;
typedef LocalOperator2D<LOp_Lapl1D_X, LOp_Id1D_Y>                                       LOp_Lapl_Id_2D;
typedef LocalOperator2D<LOp_Id1D_X, LOp_Lapl1D_Y>                                       LOp_Id_Lapl_2D;
typedef LocalOperator2D<LOp_WLapl1D_X, LOp_Id1D_Y>                                      LOp_WLapl_Id_2D;
typedef LocalOperator2D<LOp_WId1D_X, LOp_Lapl1D_Y>                                      LOp_WId_Lapl_2D;

//==== CompoundOperators ====//
typedef FlexibleCompoundLocalOperator<Index2D,AbstractLocalOperator2D<T> > 		        Flex_COp_2D;
typedef CompoundLocalOperator<Index2D,LOp_Id_Id_2D,LOp_Lapl_Id_2D,LOp_Id_Lapl_2D>       H1_InnProd_2D;

//==== Preconditioners ====//
typedef H1NormPreconditioner2D<T, Basis2D>                          Prec2D;
typedef NoPreconditioner<T, Index2D>								NoPrec2D;

//==== RightHandSides ====//
typedef SeparableRHS2D<T,Basis2D>                                   SeparableRhsIntegral2D;
typedef RHS<T,Index2D,SeparableRhsIntegral2D,
            NoPrec2D>                                         		SeparableRhs;


//==== RB Stuff ====//
const size_t PDim = 2;

typedef CoeffVector 								 						DataType;
typedef array<T,PDim>	 													ParamType;

typedef AffineLocalOperator<Index2D,AbstractLocalOperator2D<T>,ParamType>	Affine_Op_2D;
typedef AffineRhs<T,Index2D,SeparableRhs,ParamType>							Affine_Rhs_2D;
typedef FlexibleCompoundRhs<T,Index2D,SeparableRhs>							RieszF_Rhs_2D;
typedef FlexibleBilformRhs<Index2D,AbstractLocalOperator2D<T> >				RieszA_Rhs_2D;

typedef FlexibleCompoundRhs<T,Index2D,SeparableRhs>     					Flex_Rhs_2D;

typedef MultiTreeAWGM2<Index2D,Basis2D, Affine_Op_2D,
		Affine_Rhs_2D,Prec2D>                           					MT_AWGM_Truth;
typedef MultiTreeAWGM2<Index2D,Basis2D, H1_InnProd_2D,
		RieszF_Rhs_2D,Prec2D>											    MT_AWGM_Riesz_F;
typedef MultiTreeAWGM2<Index2D,Basis2D, H1_InnProd_2D,
		RieszA_Rhs_2D,Prec2D>											    MT_AWGM_Riesz_A;

typedef MT_Truth<DataType,ParamType,MT_AWGM_Truth,
				 MT_AWGM_Riesz_F,MT_AWGM_Riesz_A,H1_InnProd_2D,
				 Flex_COp_2D, Flex_Rhs_2D>									MTTruthSolver;

typedef LB_Base<ParamType, MTTruthSolver> 									LB_Type;
typedef Simple_RB_System<T,ParamType, LB_Type>								RB_Model;
typedef RB_Base<RB_Model,MTTruthSolver,DataType,ParamType>					RB_BaseModel;

#define x1 1./3.
#define x2 2./3.
#define y1 2./5.
#define y2 4./5.
#define nb_stempel 9.

T f_1_x(T x)
{
    if(x <= x1){
        return 1.;
    }
    return 0.;
}

T f_2_x(T x)
{
    if(x >= x1 && x <= x2){
        return 1.;
    }
    return 0.;
}

T f_3_x(T x)
{
    if(x >= x2){
        return 1.;
    }
    return 0.;
}

T f_1_y(T y)
{
    if(y <= y1){
        return 1.;
    }
    return 0.;
}

T f_2_y(T y)
{
    if(y >= y1 && y <= y2){
        return 1.;
    }
    return 0.;
}

T f_3_y(T y)
{
    if(y >= y2){
        return 1.;
    }
    return 0.;
}


T dummy(T, T)
{
    return 0;
}

T no_theta(const std::array<T,PDim>& /*mu*/)
{
	return 1.;
}

T theta_1(const std::array<T,PDim>& mu)
{
	return mu[0];
}

T theta_chi_1(const std::array<T,PDim>& mu)
{
    if(mu[1] < 1./nb_stempel){
       return 1;
    }
	return 0;
}

T theta_chi_2(const std::array<T,PDim>& mu)
{
    if(mu[1] >= 1./nb_stempel && mu[1] < 2./nb_stempel){
       return 1;
    }
	return 0;
}

T theta_chi_3(const std::array<T,PDim>& mu)
{
    if(mu[1] >= 2./nb_stempel && mu[1] < 3./nb_stempel){
       return 1;
    }
	return 0;
}

T theta_chi_4(const std::array<T,PDim>& mu)
{
    if(mu[1] >= 3./nb_stempel && mu[1] < 4./nb_stempel){
       return 1;
    }
	return 0;
}

T theta_chi_5(const std::array<T,PDim>& mu)
{
    if(mu[1] >= 4./nb_stempel && mu[1] < 5./nb_stempel){
       return 1;
    }
	return 0;
}

T theta_chi_6(const std::array<T,PDim>& mu)
{
    if(mu[1] >= 5./nb_stempel && mu[1] < 6./nb_stempel){
       return 1;
    }
	return 0;
}

T theta_chi_7(const std::array<T,PDim>& mu)
{
    if(mu[1] >= 6./nb_stempel && mu[1] < 7./nb_stempel){
       return 1;
    }
	return 0;
}

T theta_chi_8(const std::array<T,PDim>& mu)
{
    if(mu[1] >= 7./nb_stempel && mu[1] < 8./nb_stempel){
       return 1;
    }
	return 0;
}

T theta_chi_9(const std::array<T,PDim>& mu)
{
    if(mu[1] >= 8./nb_stempel){
       return 1;
    }
	return 0;
}

T zero_fct(T /*x*/){
	return 0;
}

T chi_omega_1(T x){
    if(x <= 0.5){
        return 1;
    }
	return 0;
}

T chi_omega_0(T x){
    if(x >= 0.5){
        return 1;
    }
	return 0;
}
