/* POISSON PROBLEM 1D
 *
 *  This example calculates a poisson problem with constant forcing f on the
 *  one-dimensional domain [0,1], i.e.
 *          - u'' = f on (0,1) , u(0) = u(1) = 0.
 *  The solution is obtained using a uniform Wavelet-Galerkin method with a
 *  diagonal scaling preconditioner.
 */

/// First we simply include the general LAWA header `lawa/lawa.h` for simplicity, thus having
/// all LAWA features available.
/// All LAWA features reside in the namespace lawa, so we introduce the `namespace lawa` globally.
#include <iostream>
#include <fstream>
#include <lawa/lawa.h>

using namespace std;
using namespace lawa;

/// Several typedefs for notational convenience.

///  Typedefs for Flens data types:
typedef double T;
typedef flens::DenseVector<flens::Array<long double> >              DenseVectorLD;

///  Typedefs for problem components:
///     Multiwavelet basis over an interval
typedef Basis<long double, Orthogonal, Interval, Multi>                       MultiWaveletBasis;
typedef Basis<long double, Orthogonal, Interval, Multi>::RefinementBasis      MultiRefinementBasis;

typedef Integral<Gauss,MultiWaveletBasis,MultiWaveletBasis>         MultiWaveletIntegral;
typedef Integral<Gauss,MultiRefinementBasis,MultiRefinementBasis>   MultiRefinementIntegral;

int main(int argc, char*argv[])
{
    cout.precision(20);
    /// wavelet basis parameters:
    int d = 4;
    int j0 = 3;
    if (d==3) { j0 = 3; }
    int J = 3;

    /// Basis initialization, using Dirichlet boundary conditions
    MultiWaveletBasis multibasis(d, 0);     // For L2_orthonormal and special MW bases
    multibasis.enforceBoundaryCondition<DirichletBC>();

/*
    cout << 1./std::sqrt(2.) << " " << pow2ih<long double>(-1) << " " << endl;
    cout << 2.L/3.L << " " << _quadratic_refinement_left_evaluator0<long double>(2.L/3.L,0) << endl;
    cout << 2.L/3.L << " " << pow2ih<long double>(-1)*multibasis.refinementbasis.generator(XBSpline).operator()(1.L/3.L,1,0,0) << endl;
    cout << 2.L/3.L << " " << 1.L/std::sqrt(2.L)*multibasis.refinementbasis.generator(XBSpline).operator()(1.L/3.L,1,0,0) << endl << endl;

    cout << 2.L/3.L << " " << _quadratic_refinement_right_evaluator0<long double>(4.L/3.L,0) << endl;
*/
    /*
    for (int j=1; j<=J; ++j) {
        cout << "j = " << j << ": " << multibasis.refinementbasis.mra.cardI(j) << " " << multibasis.refinementbasis.mra.rangeI(j) << endl;
        for (int k= multibasis.refinementbasis.mra.rangeI(j).firstIndex();
                 k<=multibasis.refinementbasis.mra.rangeI(j).lastIndex(); ++k) {
            cout << "  k = " << k  << " "
                 << multibasis.refinementbasis.generator(XBSpline).support(j,k) << " "
                 << multibasis.refinementbasis.generator(XBSpline).singularSupport(j,k) << endl;
            ofstream plotfile_scaling("refinement_interval.txt");
            for (T x=0.; x<=1.; x+=pow2i<T>(-6-j)) {
                plotfile_scaling << x << " " << multibasis.refinementbasis.generator(XBSpline).operator()(x,j,k,0) << endl;
            }
            plotfile_scaling.close();
            getchar();
        }
    }
    */

    /// Test refinement of inner multi scaling functions

    DenseVectorLD *refCoeffs;
    int deriv = 1;
/*
    for (int j=0; j<=J; ++j) {
        for (int k=multibasis.mra.rangeI(j).firstIndex(); k<=multibasis.mra.rangeI(j).lastIndex(); ++k) {
            ofstream plotfile_scaling("refinement_interval_multiscaling.txt");

            int refinement_j = 0;
            long refinement_k_first = 0L;
            refCoeffs = multibasis.mra.phi.getRefinement(j,k,refinement_j,refinement_k_first);
            cout << "j = " << j << ", k = " << k << ", refinement_j = "
                           << refinement_j << ", refinement_k_first = " << refinement_k_first << endl;
            for (T x=0.; x<=1.; x+=pow2i<T>(-12-j)) {
                T reference_value = multibasis.generator(XBSpline)(x,j,k,0);
                T refinement_value = 0.;
                for (int i=(*refCoeffs).firstIndex(); i<=(*refCoeffs).lastIndex(); ++i) {
                    refinement_value +=   (*refCoeffs).operator()(i)
                                        * multibasis.refinementbasis.generator(XBSpline).operator()(x,refinement_j,refinement_k_first+i,0);
                }
                plotfile_scaling << x << " " << reference_value << " " << refinement_value << endl;
            }
            plotfile_scaling.close();
            cout << "Please hit enter." << endl;
            getchar();
        }
    }
*/

    for (int j=0; j<=J; ++j) {
        for (int k=multibasis.rangeJ(j).firstIndex(); k<=multibasis.rangeJ(j).lastIndex(); ++k) {
            ofstream plotfile_wavelet("refinement_interval_multiwavelet.txt");
            long double error = 0.L;
            int refinement_j = 0;
            long refinement_k_first = 0L;
            refCoeffs = multibasis.psi.getRefinement(j,k,refinement_j,refinement_k_first);
            cout << "j = " << j << ", k = " << k << ", refinement_j = " << refinement_j
                 << ", refinement_k_first = " << refinement_k_first << endl;
            for (T x=0.; x<=1.; x+=pow2i<T>(-12-j)) {
                long double reference_value = multibasis.generator(XWavelet)(x,j,k,deriv);
                long double refinement_value = 0.L;
                for (int i=(*refCoeffs).firstIndex(); i<=(*refCoeffs).lastIndex(); ++i) {
                    refinement_value +=   (long double)(*refCoeffs).operator()(i)
                                        * multibasis.refinementbasis.generator(XBSpline).operator()(x,refinement_j,refinement_k_first+i,deriv);
                }
                error = std::max(error, fabs(reference_value-refinement_value)/*/fabs(reference_value)*/);
                plotfile_wavelet << x << " " << reference_value << " " << refinement_value << endl;
            }
            plotfile_wavelet.close();
            cout << "Error in || ||_infty = " << error << endl;
            cout << "Please hit enter." << endl;
            getchar();
        }
    }



    MultiRefinementIntegral   refinement_integral(multibasis.refinementbasis,multibasis.refinementbasis);
    MultiWaveletIntegral multi_integral(multibasis,multibasis);
    int j=3;
    DenseVectorLD *refCoeffs1, *refCoeffs2;
    int refinement_j1 = 0, refinement_j2 = 0;
    long refinement_k_first1 = 0L, refinement_k_first2 = 0L;
    refCoeffs1 = multibasis.mra.phi.getRefinement(j,2,refinement_j1,refinement_k_first1);
    refCoeffs2 = multibasis.mra.phi.getRefinement(j,3,refinement_j2,refinement_k_first2);

    long double val = 0.;
    for (int i1=(*refCoeffs1).firstIndex(); i1<=(*refCoeffs1).lastIndex(); ++i1) {
        for (int i2=(*refCoeffs2).firstIndex(); i2<=(*refCoeffs2).lastIndex(); ++i2) {
            T integral_value = refinement_integral(refinement_j1,refinement_k_first1+i1,XBSpline,0,
                                                   refinement_j2,refinement_k_first2+i2,XBSpline,0);
            cout << (*refCoeffs1).operator()(i1)*(*refCoeffs2).operator()(i2)*integral_value << endl;
            val += (*refCoeffs1).operator()(i1)*(*refCoeffs2).operator()(i2)*integral_value;
        }
    }
    cout << "Value = " << val << " " << multi_integral(j,2,XBSpline,0, j,3,XBSpline,0) << endl;

    return 0;
}
