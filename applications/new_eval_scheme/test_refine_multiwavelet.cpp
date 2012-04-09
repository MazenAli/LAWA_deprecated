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
typedef long double T;
typedef flens::DenseVector<flens::Array<long double> >              DenseVectorLD;

///  Typedefs for problem components:
///     Multiwavelet basis over an interval
typedef Basis<T, Orthogonal, Interval, Multi>                       MultiWaveletBasis;
typedef Basis<T, Orthogonal, Interval, Multi>::RefinementBasis      MultiRefinementBasis;

typedef Integral<Gauss,MultiWaveletBasis,MultiWaveletBasis>         MultiWaveletIntegral;
typedef Integral<Gauss,MultiRefinementBasis,MultiRefinementBasis>   MultiRefinementIntegral;

int main(int argc, char*argv[])
{
    cout.precision(20);
    /// wavelet basis parameters:
    int d = 3;
    int j0 = 3;
    if (d==3) { j0 = 3; }
    int J = 3;

    /// Basis initialization, using Dirichlet boundary conditions
    MultiWaveletBasis multibasis(d, 0);     // For L2_orthonormal and special MW bases
    multibasis.enforceBoundaryCondition<DirichletBC>();

    /// Check precision of some generators
    cout << 1./std::sqrt(2.L) << " " << pow2ih<long double>(-1) << " " << endl;
    cout << 2.L/3.L << " " << _quadratic_refinement_left_evaluator0<long double>(2.L/3.L,0) << endl;
    cout << 2.L/9.L << " " << _quadratic_refinement_left_evaluator0<long double>(4.L/3.L,0) << endl << endl;
    cout << 4.L/9.L << " " << _cubic_refinement_left_evaluator0<long double>(1.L/3.L,0) << endl;
    cout << 2.L/9.L << " " << _cubic_refinement_right_evaluator0<long double>(4.L/3.L,0) << endl << endl;

    /*
    /// Plot refinement basis
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
            T abs_error = 0.L, rel_error = 0.L;
            T x_rel_crit = 0.L, x_abs_crit = 0.L;
            int refinement_j = 0;
            long refinement_k_first = 0L;
            refCoeffs = multibasis.mra.phi.getRefinement(j,k,refinement_j,refinement_k_first);
            cout << "j = " << j << ", k = " << k << ", refinement_j = "
                           << refinement_j << ", refinement_k_first = " << refinement_k_first << endl;
            for (T x=0.; x<=1.; x+=pow2i<T>(-8-j)) {
                T reference_value = multibasis.generator(XBSpline)(x,j,k,deriv);
                T refinement_value = 0.;
                for (int i=(*refCoeffs).firstIndex(); i<=(*refCoeffs).lastIndex(); ++i) {
                    refinement_value +=   (T)(*refCoeffs).operator()(i)
                                        * multibasis.refinementbasis.generator(XBSpline).operator()(x,refinement_j,refinement_k_first+i,deriv);
                }
                if (x==0.5) {
                    cout << "x = 0.5: " << reference_value << " " << refinement_value << endl;
                }
                T tmp1 = fabs(reference_value-refinement_value)/fabs(reference_value);
                T tmp2 = fabs(reference_value-refinement_value);
                if (rel_error<=tmp1 && fabs(reference_value)>1e-15) { rel_error = tmp1; x_rel_crit = x;    }
                if (abs_error<=tmp2) { abs_error = tmp2; x_abs_crit = x;    }
                plotfile_scaling << x << " " << reference_value << " " << refinement_value << endl;
            }
            plotfile_scaling.close();
            cout << "Relative error: " << rel_error << ", x_rel_crit = " << x_rel_crit << endl;
            cout << "Absolute error: " << abs_error << ", x_abs_crit = " << x_abs_crit << endl;
            cout << "Please hit enter." << endl;
            getchar();
        }
    }

    for (int j=0; j<=J; ++j) {
        for (int k=multibasis.rangeJ(j).firstIndex(); k<=multibasis.rangeJ(j).lastIndex(); ++k) {
            ofstream plotfile_wavelet("refinement_interval_multiwavelet.txt");
            T abs_error = 0.L, rel_error = 0.L;
            T x_rel_crit = 0.L, x_abs_crit = 0.L;
            int refinement_j = 0;
            long refinement_k_first = 0L;
            refCoeffs = multibasis.psi.getRefinement(j,k,refinement_j,refinement_k_first);
            cout << "j = " << j << ", k = " << k << ", refinement_j = " << refinement_j
                 << ", refinement_k_first = " << refinement_k_first << endl;
            for (T x=0.L; x<=1.L; x+=pow2i<long double>(-8-j)) {
                T reference_value = multibasis.generator(XWavelet)(x,j,k,deriv);
                T refinement_value = 0.L;
                for (int i=(*refCoeffs).firstIndex(); i<=(*refCoeffs).lastIndex(); ++i) {
                    refinement_value +=   (T)(*refCoeffs).operator()(i)
                                        * multibasis.refinementbasis.generator(XBSpline).operator()(x,refinement_j,refinement_k_first+i,deriv);
                }
                T tmp1 = fabs(reference_value-refinement_value)/fabs(reference_value);
                T tmp2 = fabs(reference_value-refinement_value);
                if (rel_error<=tmp1 && fabs(reference_value)>1e-15) { rel_error = tmp1; x_rel_crit = x;    }
                if (abs_error<=tmp2) { abs_error = tmp2; x_abs_crit = x;    }
                plotfile_wavelet << x << " " << reference_value << " " << refinement_value << endl;
            }
            plotfile_wavelet.close();
            cout << "Relative error: " << rel_error << ", x_rel_crit = " << x_rel_crit << endl;
            cout << "Absolute error: " << abs_error << ", x_abs_crit = " << x_abs_crit << endl;
            cout << "Please hit enter." << endl;
            getchar();
        }
    }
    */

    MultiRefinementIntegral   refinement_integral(multibasis.refinementbasis,multibasis.refinementbasis);
    MultiWaveletIntegral multi_integral(multibasis,multibasis);
    DenseVectorLD *refCoeffs1, *refCoeffs2;
    T max_error = 0.L;
    for (int j1=0; j1<=J; ++j1) {
    for (int k1=multibasis.rangeJ(j1).firstIndex(); k1<=multibasis.rangeJ(j1).lastIndex(); ++k1) {
        for (int j2=0; j2<=J; ++j2) {
        for (int k2=multibasis.rangeJ(j2).firstIndex(); k2<=multibasis.rangeJ(j2).lastIndex(); ++k2) {
            int refinement_j1 = 0, refinement_j2 = 0;
            long refinement_k_first1 = 0L, refinement_k_first2 = 0L;
            refCoeffs1 = multibasis.psi.getRefinement(j1,k1,refinement_j1,refinement_k_first1);
            refCoeffs2 = multibasis.psi.getRefinement(j2,k2,refinement_j2,refinement_k_first2);
            T val = 0.L;
            for (int i1=(*refCoeffs1).firstIndex(); i1<=(*refCoeffs1).lastIndex(); ++i1) {
                for (int i2=(*refCoeffs2).firstIndex(); i2<=(*refCoeffs2).lastIndex(); ++i2) {
                    T integral_value = refinement_integral(refinement_j1,refinement_k_first1+i1,XBSpline,0,
                                                           refinement_j2,refinement_k_first2+i2,XBSpline,0);
                    val += (*refCoeffs1).operator()(i1)*(*refCoeffs2).operator()(i2)*integral_value;
                }
            }
            if (j1==j2 && k1==k2) {
                cout << "(" << j1 << "," << k1 << "), (" << j2 << "," << k2 << "), Error: " << fabs(val-1.) << endl;
                max_error = std::max(max_error,fabs(val-1.));
            }
            else  {
                cout << "(" << j1 << "," << k1 << "), (" << j2 << "," << k2 << "), Error: " << fabs(val) << endl;
                max_error = std::max(max_error,fabs(val));
            }

        }
        }
    }
    }
    cout << "Max integration error: " << max_error << endl;

    return 0;
}
