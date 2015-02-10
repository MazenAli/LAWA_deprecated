//
//  problem.h
//  LAWA
//
//  Created by Silke Glas.
//  04.02.2012 Erweiterung der problem.h aus ex1 um output
//

#include <iostream>
#include <iomanip>
#include <lawa/lawa.h>
#include <lawa/methods/rb/postprocessing/plotting.h>
#include <lawa/methods/rb/deprecated/datastructures/datastructures.h>
#include <lawa/methods/rb/deprecated/outputs/averageoutput.h>
#include <lawa/methods/rb/deprecated/solvers/solvers.h>
/*********************************************************************
 *    EXAMPLE
 *********************************************************************/
 
/* Example: Thermal block  
 *
 *          - theta * u_xx = f on (0,1)^2
 *         u(x,0) = u(x,1) = 0
 *
 *        Here: theta = theta_1 on (0, 0.5) x (0,1)
 *              theta = theta_2 on (0,5, 1) x (0,1)
 *                  f = 1       on (0.4, 0.6) x (0.45, 0.55)
 */
 
/*********************************************************************
 *    TYPEDEFS
 *********************************************************************/
 

using namespace std;
using namespace lawa;

typedef double T;

// Basis definitions
typedef Basis<T,Primal,Interval,Dijkema>                    IntervalBasis;
typedef TensorBasis2D<Adaptive,IntervalBasis,IntervalBasis> Basis2D;

// Operator Definitions
typedef WeightedHelmholtzOperator2D<T, Basis2D>                           WeightedLaplaceOp2D;
typedef H1NormPreconditioner2D<T, Basis2D>                                DiagPrec2D;
typedef NoPreconditioner<T, Index2D>                                    NoPrec2D;
typedef WeightedAdaptiveHelmholtzOperator2D<T, Basis2D, NoPrec2D>         WeightedAdaptHHOp2D;
typedef AdaptiveHelmholtzOperator2D<T, Basis2D, NoPrec2D>                 AdaptHHOp2D;
typedef RHS<T,Index2D, SeparableRHS2D<T, Basis2D>, NoPrec2D>              AdaptRHS;
typedef AverageOutput2D<T, Index2D, Basis2D>							 Output;

// Algorithm Definition
// Class for calling the truth solver for snapshot calculations
typedef CompressionWeightedPDE2D<T, Basis2D>                                    Compression;
typedef IndexsetTruthSolver<T, Basis2D, DiagPrec2D, Index2D, Compression>       IndexsetSolver;
// Class containing all \calN-dependent data and functions
typedef AdaptiveRBTruth2D<T, Basis2D, DiagPrec2D,IndexsetSolver, Compression>   RBTruth;
// Class containing only N-dependent data and functions
typedef RBModel2D<T, RBTruth>                                                   RBModel;

// Data type definitions
typedef Coefficients<Lexicographical,T,Index2D>                       Coeffs;
typedef flens::DenseVector<flens::Array<T> >                          DenseVectorT;
typedef flens::GeMatrix<flens::FullStorage<T, cxxblas::ColMajor> >    FullColMatrixT;


/*********************************************************************
 *    PARAMETERS AND CONSTANTS
 *********************************************************************/

const double reference_mu = 1;
const double mu_min = 0.1;
const double mu_max = 2.;
const double mu_init = 2.;

const int    n_train = 30;
const double Greedy_tol = 1e-10;
const int    Nmax = 7;

const double test_min = 0.13;
const double test_max = 0.96;
const int    n_test = 29;
const double test_h = (mu_max - mu_min)/(n_test-1);

/*********************************************************************
 *    AFFINE DECOMPOSITIONS
 *********************************************************************/

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
weight_Omega_output(T output)
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
theta_f_1(const std::vector<T>&)
{
    return 1.;
}

T
theta_output_1(const std::vector<T>&)
{
    return 1.;
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
plot_dummy_fct(T x, T y)
{
    return 0;
}
