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
typedef Basis<T, Orthogonal, Interval, Multi>                       MultiWaveletBasis;
typedef Basis<T, Orthogonal, Interval, Multi>::RefinementBasis      MultiRefinementBasis;

typedef Integral<Gauss,MultiWaveletBasis,MultiWaveletBasis>         MultiWaveletIntegral;
typedef Integral<Gauss,MultiRefinementBasis,MultiRefinementBasis>   MultiRefinementIntegral;

int main(int argc, char*argv[])
{
    /// wavelet basis parameters:
    int d = 3;
    int j0 = 3;
    if (d==3) { j0 = 3; }
    int J = 3;

    /// Basis initialization, using Dirichlet boundary conditions
    MultiWaveletBasis multibasis(d, 0);     // For L2_orthonormal and special MW bases
    multibasis.enforceBoundaryCondition<DirichletBC>();

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


    /// Test refinement of inner multi scaling functions
    /*
    DenseVectorLD *refCoeffs;
    int addRefLevel = multibasis.mra.addRefLevel;
    int factor   = pow2i<int>(addRefLevel);


    for (int j=0; j<=J; ++j) {
        for (int k=multibasis.mra.rangeI(j).firstIndex(); k<=multibasis.mra.rangeI(j).lastIndex(); ++k) {
            ofstream plotfile_scaling("refinement_interval_multiscaling.txt");

            long shift = 0L;
            long offset = 0L;
            refCoeffs = multibasis.mra.phi.getRefinement(j,k,shift,offset);
            cout << "j = " << j << ", k = " << k << ", shift = " << shift << ", offset = " << offset << endl;
            for (T x=0.; x<=1.; x+=pow2i<T>(-6-j)) {
                T reference_value = multibasis.generator(XBSpline)(x,j,k,0);
                T refinement_value = 0.;
                for (int i=(*refCoeffs).firstIndex(); i<=(*refCoeffs).lastIndex(); ++i) {
                    refinement_value +=   (*refCoeffs).operator()(i)
                                        * multibasis.refinementbasis.generator(XBSpline).operator()(x,j+addRefLevel,i+factor*shift+offset,0);
                }
                plotfile_scaling << x << " " << reference_value << " " << refinement_value << endl;
            }
            plotfile_scaling.close();
            cout << "Please hit enter." << endl;
            getchar();
        }
    }
    for (int j=0; j<=J; ++j) {
        for (int k=multibasis.rangeJ(j).firstIndex(); k<=multibasis.rangeJ(j).lastIndex(); ++k) {
            ofstream plotfile_wavelet("refinement_interval_multiwavelet.txt");

            long shift = 0L;
            long offset = 0L;
            refCoeffs = multibasis.psi.getRefinement(j,k,shift,offset);
            cout << "j = " << j << ", k = " << k << ", shift = " << shift << ", offset = " << offset << endl;
            for (T x=0.; x<=1.; x+=pow2i<T>(-6-j)) {
                T reference_value = multibasis.generator(XWavelet)(x,j,k,0);
                T refinement_value = 0.;
                for (int i=(*refCoeffs).firstIndex(); i<=(*refCoeffs).lastIndex(); ++i) {
                    refinement_value +=   (*refCoeffs).operator()(i)
                                        * multibasis.refinementbasis.generator(XBSpline).operator()(x,j+addRefLevel,i+factor*shift+offset,0);
                }
                plotfile_wavelet << x << " " << reference_value << " " << refinement_value << endl;
            }
            plotfile_wavelet.close();
            cout << "Please hit enter." << endl;
            getchar();
        }
    }

    */
    /*
    RefinementIntegral   refinement_integral(refinebasis,refinebasis);
    MultiWaveletIntegral multi_integral(multibasis,multibasis);
    int j=3;
    DenseVectorLD *refCoeffs1, *refCoeffs2;
    long shift1 = 0L, offset1 = 0L, shift2 = 0L, offset2 = 0L;
    refCoeffs1 = multibasis.mra.phi.getRefinement(j,2,shift1,offset1);
    refCoeffs2 = multibasis.mra.phi.getRefinement(j,3,shift2,offset2);

    long double val = 0.;
    for (int i1=(*refCoeffs1).firstIndex(); i1<=(*refCoeffs1).lastIndex(); ++i1) {
        int k1 = i1+factor*shift1+offset1;
        for (int i2=(*refCoeffs2).firstIndex(); i2<=(*refCoeffs2).lastIndex(); ++i2) {
            int k2 = i2+factor*shift2+offset2;
            T integral_value = refinement_integral(j+addRefLevel,k1,XBSpline,0, j+addRefLevel,k2,XBSpline,0);
            cout << (*refCoeffs1).operator()(i1)*(*refCoeffs2).operator()(i2)*integral_value << endl;
            val += (*refCoeffs1).operator()(i1)*(*refCoeffs2).operator()(i2)*integral_value;
            cout << "  --> val = " << val << endl;
        }
    }
    cout << "Value = " << val << " " << multi_integral(j,2,XBSpline,0, j,3,XBSpline,0) << endl;
    */
    return 0;
}
