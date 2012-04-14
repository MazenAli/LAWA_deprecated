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
//typedef Basis<T, Orthogonal, Interval, Multi>                       PrimalBasis;
typedef Basis<T, Primal, Interval, Dijkema>                       PrimalBasis;
typedef PrimalBasis::RefinementBasis                                RefinementBasis;

typedef Integral<Gauss,PrimalBasis,PrimalBasis>                     MultiWaveletIntegral;
typedef Integral<Gauss,RefinementBasis,RefinementBasis>             MultiRefinementIntegral;

int main(int argc, char*argv[])
{
    cout.precision(20);
    if (argc!=4) {
        cout << "Usage: " << argv[0] << " d j0 J" << endl;
        return 0;
    }
    /// wavelet basis parameters:
    int d  = atoi(argv[1]);
    int j0 = atoi(argv[2]);
    int J  = atoi(argv[3]);
    int deriv = 0;

    /// refinement coefficient vector
    DenseVectorLD *refCoeffs;


    /// Basis initialization, using Dirichlet boundary conditions
    //PrimalBasis basis(d, j0);     // For L2_orthonormal and special MW bases
    PrimalBasis basis(d, d, j0);     // For biorthogonal wavelet bases
    basis.enforceBoundaryCondition<DirichletBC>();
    RefinementBasis& refinementbasis = basis.refinementbasis;

    /// Test refinement of B-splines
    cout << " *******************************************" << endl;
    cout << " ******** Refinement of B-splines **********" << endl;
    for (int j=refinementbasis.mra.j0; j<=refinementbasis.mra.j0+2; ++j) {
        cout << "j = " << j << ": " << refinementbasis.mra.cardI(j) << " " << refinementbasis.mra.rangeI(j) << endl;
        for (int k= refinementbasis.mra.rangeI(j).firstIndex(); k<=refinementbasis.mra.rangeI(j).lastIndex(); ++k) {
            cout << "  k = " << k  << " " << refinementbasis.generator(XBSpline).support(j,k) << " "
                 << refinementbasis.generator(XBSpline).singularSupport(j,k) << endl;
            ofstream plotfile_scaling("refinement_interval.txt");
            T abs_error = 0.L, rel_error = 0.L;
            T x_rel_crit = 0.L, x_abs_crit = 0.L;
            int refinement_j = 0;
            long refinement_k_first = 0L;
            refCoeffs = refinementbasis.mra.phi.getRefinement(j,k,refinement_j,refinement_k_first);
            for (T x=0.; x<=1.; x+=pow2i<T>(-10-j)) {
                T reference_value = refinementbasis.generator(XBSpline).operator()(x,j,k,deriv);
                T refinement_value = 0.L;
                for (int i=(*refCoeffs).firstIndex(); i<=(*refCoeffs).lastIndex(); ++i) {
                    refinement_value +=   (T)(*refCoeffs).operator()(i)
                                        * refinementbasis.generator(XBSpline).operator()(x,refinement_j,refinement_k_first+i,deriv);
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
    cout << " *******************************************" << endl << endl;

    /// Test refinement of multiscaling functions
    cout << " *******************************************" << endl;
    cout << " ******* Refinement of multiscalings *******" << endl;
    for (int j=j0; j<=J; ++j) {
        for (int k=basis.mra.rangeI(j).firstIndex(); k<=basis.mra.rangeI(j).lastIndex(); ++k) {
            ofstream plotfile_scaling("refinement_interval_multiscaling.txt");
            T abs_error = 0.L, rel_error = 0.L;
            T x_rel_crit = 0.L, x_abs_crit = 0.L;
            int refinement_j = 0;
            long refinement_k_first = 0L;
            refCoeffs = basis.mra.phi.getRefinement(j,k,refinement_j,refinement_k_first);
            cout << "j = " << j << ", k = " << k << ", refinement_j = "
                           << refinement_j << ", refinement_k_first = " << refinement_k_first << endl;
            cout << refinementbasis.mra.rangeI(refinement_j) << endl;
            for (T x=0.; x<=1.; x+=pow2i<T>(-8-j)) {
                T reference_value = basis.generator(XBSpline)(x,j,k,deriv);
                T refinement_value = 0.;
                for (int i=(*refCoeffs).firstIndex(); i<=(*refCoeffs).lastIndex(); ++i) {
                    refinement_value +=   (T)(*refCoeffs).operator()(i)
                                        * basis.refinementbasis.generator(XBSpline).operator()(x,refinement_j,refinement_k_first+i,deriv);
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
    cout << " *******************************************" << endl << endl;

    /// Test refinement of multiwavelets
    cout << " *******************************************" << endl;
    cout << " ******* Refinement of multiwavelets *******" << endl;
    for (int j=j0; j<=J; ++j) {
        for (int k=basis.rangeJ(j).firstIndex(); k<=basis.rangeJ(j).lastIndex(); ++k) {
            ofstream plotfile_wavelet("refinement_interval_multiwavelet.txt");
            T abs_error = 0.L, rel_error = 0.L;
            T x_rel_crit = 0.L, x_abs_crit = 0.L;
            int refinement_j = 0;
            long refinement_k_first = 0L;
            refCoeffs = basis.psi.getRefinement(j,k,refinement_j,refinement_k_first);
            long refinement_k_last  = 0L;
            basis.psi.getRefinementNeighbors(j,k,refinement_j,refinement_k_first,refinement_k_last);
            cout << "j = " << j << ", k = " << k << ", refinement_j = " << refinement_j
                 << ", refinement_k_first = " << refinement_k_first << endl;
            cout << "Refinement range: [" << refinement_k_first << ", " << refinement_k_first+(*refCoeffs).lastIndex() << "]" << endl;
            cout << "                  [" << refinement_k_first << ", " << refinement_k_last << "]" << endl;
            for (T x=0.L; x<=1.L; x+=pow2i<long double>(-8-j)) {
                T reference_value = basis.generator(XWavelet)(x,j,k,deriv);
                T refinement_value = 0.L;
                for (int i=(*refCoeffs).firstIndex(); i<=(*refCoeffs).lastIndex(); ++i) {

                    refinement_value +=   (T)(*refCoeffs).operator()(i)
                                        * refinementbasis.generator(XBSpline).operator()(x,refinement_j,refinement_k_first+i,deriv);
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


    /// Check for refined neighbors: Given a b-spline, we need to determine the wavelets whose
    /// supports intersect the one of the b-spline. Concerning the levels: suppose that j indicates
    /// the level of $\psi_{j,k}$. Then we reconstruct $\psi_{j,k}$ by b-splines on j_refinement+1.
    /// Note: the level j_refinement corresponds to the level of b-splines required for the refinement
    /// of a multiscaling function $\phi_{j,k}$.
    for (int refinement_j=j0+4; refinement_j<=J+4; ++refinement_j) {
        for (int refinement_k =refinementbasis.mra.rangeI(refinement_j).firstIndex();
                 refinement_k<=refinementbasis.mra.rangeI(refinement_j).lastIndex(); ++refinement_k) {
            int j=0;
            long k_first=0L, k_last=0L;
            cout << "(" << refinement_j << "," << refinement_k << "): " << endl;
            basis.psi.getRefinedNeighbors(refinement_j, refinement_k, j, k_first, k_last);
            cout << "(" << refinement_j << "," << refinement_k << "): "
                             << j << " , [" << k_first << "," << k_last << "], "
                             << basis.rangeJ(j) << endl;;

            for (long k=basis.rangeJ(j).firstIndex(); k<k_first; ++k) {
                if (overlap(refinementbasis.mra.phi.support(refinement_j,refinement_k),
                            basis.psi.support(j,k))>0) {
                    cout << "Error: k=" << k << " in " << basis.rangeJ(j) << " is missing."
                         << refinementbasis.mra.phi.support(refinement_j,refinement_k)
                         << " " << basis.psi.support(j,k) << endl;
                }
            }
            for (long k=k_last+1; k<basis.rangeJ(j).lastIndex(); ++k) {
                if (overlap(refinementbasis.mra.phi.support(refinement_j,refinement_k),
                            basis.psi.support(j,k))>0) {
                    cout << "Error: k=" << k << " in " << basis.rangeJ(j) << " is missing."
                         << refinementbasis.mra.phi.support(refinement_j,refinement_k)
                         << " " << basis.psi.support(j,k) << endl;
                }
            }
            cout << endl;
            getchar();
        }
    }

    /*
    MultiRefinementIntegral   refinement_integral(basis.refinementbasis,basis.refinementbasis);
    MultiWaveletIntegral multi_integral(basis,basis);
    DenseVectorLD *refCoeffs1, *refCoeffs2;
    T max_error = 0.L;
    for (int j1=0; j1<=J; ++j1) {
    for (int k1=basis.rangeJ(j1).firstIndex(); k1<=basis.rangeJ(j1).lastIndex(); ++k1) {
        for (int j2=0; j2<=J; ++j2) {
        for (int k2=basis.rangeJ(j2).firstIndex(); k2<=basis.rangeJ(j2).lastIndex(); ++k2) {
            int refinement_j1 = 0, refinement_j2 = 0;
            long refinement_k_first1 = 0L, refinement_k_first2 = 0L;
            refCoeffs1 = basis.psi.getRefinement(j1,k1,refinement_j1,refinement_k_first1);
            refCoeffs2 = basis.psi.getRefinement(j2,k2,refinement_j2,refinement_k_first2);
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
    */
    return 0;
}


/*
    /// Check precision of some generators
    long double OneDivSqrt2 = 0.707106781186547524401L;
    cout << 1./std::sqrt(2.L)-OneDivSqrt2 << " " << std::pow(2.L,-0.5L)-OneDivSqrt2
         << " " << pow2ih<long double>(-1)-OneDivSqrt2 << " " << pow2ih<double>(-1)-OneDivSqrt2 << endl;
    cout << 2.L/3.L << " " << _quadratic_refinement_left_evaluator0<long double>(2.L/3.L,0) << endl;
    cout << 2.L/9.L << " " << _quadratic_refinement_left_evaluator0<long double>(4.L/3.L,0) << endl << endl;
    cout << 4.L/9.L << " " << _cubic_refinement_left_evaluator0<long double>(1.L/3.L,0) << endl;
    cout << 2.L/9.L << " " << _cubic_refinement_right_evaluator0<long double>(4.L/3.L,0) << endl << endl;
*/

