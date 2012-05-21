/* TEST LOCAL DECOMPOSITION
 *
 */

#include <iostream>
#include <fstream>
#include <lawa/lawa.h>

using namespace std;
using namespace lawa;

/// Several typedefs for notational convenience.

///  Typedefs for Flens data types:
typedef long double T;
typedef flens::DenseVector<flens::Array<T> >                        DenseVectorT;
typedef flens::GeMatrix<flens::FullStorage<T, cxxblas::ColMajor> >  DenseMatrixT;

///  Typedefs for problem components:
///  Wavelet basis over an interval
typedef Basis<T, Orthogonal, Interval, Multi>                       PrimalBasis;
//typedef Basis<T, Primal, Interval, Dijkema>                         PrimalBasis;
typedef PrimalBasis::RefinementBasis                                RefinementBasis;

///  Wavelet integrals
typedef IntegralF<Gauss,PrimalBasis>                                IntegralFBasis;
typedef IntegralF<Gauss,RefinementBasis>                            IntegralFRefinementBasis;

typedef CoefficientsByLevel<T>::const_it                            const_coeffbylevel_it;
typedef IndexSet<Index1D>::const_iterator                           const_set1d_it;
typedef Coefficients<Lexicographical,T,Index1D>::const_iterator     const_coeff1d_it;

T
f(T x) { return exp(x); }

int main(int argc, char*argv[])
{
    cout.precision(20);
    if (argc!=3) {
        cout << "Usage: " << argv[0] << " d j0" << endl;
        return 0;
    }
    /// wavelet basis parameters:
    int d  = atoi(argv[1]);
    int j0 = atoi(argv[2]);
    int j_wavelet  = j0+5;
    int j_scaling  = j_wavelet;

    /// Basis initialization, using Dirichlet boundary conditions
    PrimalBasis basis(d, j0);           // For L2_orthonormal and special MW bases
    //PrimalBasis basis(d, d, j0);      // For biorthogonal wavelet bases
    basis.enforceBoundaryCondition<DirichletBC>();
    RefinementBasis &refinementbasis = basis.refinementbasis;

    /// Integral initialization
    DenseVectorT                singPts;
    Function<T>                 f_fct(f,singPts);
    IntegralFBasis              integralf_basis(f_fct,basis);
    IntegralFRefinementBasis    integralf_refinementbasis(f_fct,refinementbasis);
    integralf_basis.quadrature.setOrder(30);
    integralf_refinementbasis.quadrature.setOrder(30);

    /// Local refinement initialization
    LocalRefinement<PrimalBasis> LocalRefine(basis);

    /// Computing a vector of refinement B-spline coefficients
    CoefficientsByLevel<T> f_loc_single;
    int j_refinement = basis.psi.getRefinementLevel(j_wavelet);
    for (int k= refinementbasis.mra.rangeI(j_refinement).firstIndex();
             k<=refinementbasis.mra.rangeI(j_refinement).lastIndex(); ++k) {
        f_loc_single.map.operator[](k) = integralf_refinementbasis(j_refinement,k,XBSpline,0);
    }

    /// Initializing output vectors for (multi-)wavelets and (multi-)scalings.
    CoefficientsByLevel<T> f_bspline, f_wavelet, f_scaling;
    for (int k =basis.rangeJ(j_wavelet).firstIndex(); k<=basis.rangeJ(j_wavelet).lastIndex(); ++k) {
        f_wavelet.map.operator[](k) = 0.;
    }
    for (int k =refinementbasis.mra.rangeI(j_refinement-1).firstIndex();
             k<=refinementbasis.mra.rangeI(j_refinement-1).lastIndex(); ++k) {
        f_bspline.map.operator[](k) = 0.;
    }
    if (PrimalBasis::Cons==Multi) {
        for (int k =basis.mra.rangeI(j_scaling).firstIndex();
                 k<=basis.mra.rangeI(j_scaling).lastIndex(); ++k) {
            f_scaling.map.operator[](k) = 0.;
        }
    }

    /// Decompose the refinement B-spline coefficient vector.
    LocalRefine.decompose_(f_loc_single, f_bspline, j_refinement-1, f_wavelet, j_wavelet);

    cout << "************** Wavelet Decompositions **************" << endl;
    for (const_coeffbylevel_it it=f_wavelet.map.begin(); it!=f_wavelet.map.end(); ++it) {
        T val1 = integralf_basis(j_wavelet,(*it).first,XWavelet,0);
        T val2 = (*it).second;
        cout << (*it).first << ": " << val1 << " " << val2 << " " << fabs(val1-val2) << " "
             << fabs(val1-val2)/fabs(val1) << endl;
    }
    cout << endl;
    cout << "************** B-spline Decompositions *************" << endl;
    for (const_coeffbylevel_it it=f_bspline.map.begin(); it!=f_bspline.map.end(); ++it) {
        T val1 = integralf_refinementbasis(j_refinement-1,(*it).first,XBSpline,0);
        T val2 = (*it).second;
        cout << (*it).first << ": " << val1 << " " << val2 << " " << fabs(val1-val2) << " "
             << fabs(val1-val2)/fabs(val1) << endl;
    }
    cout << endl;

    /// In case of multiscalings, we have to further decompose the refinement B-splines.
    if (PrimalBasis::Cons==Multi) {
        LocalRefine.decompose_OnlyMultiScaling(f_bspline, f_scaling, j_scaling);
        cout << "************** Scaling Decompositions *************" << endl;
        for (const_coeffbylevel_it it=f_scaling.map.begin(); it!=f_scaling.map.end(); ++it) {
            T val1 = integralf_basis(j_scaling,(*it).first,XBSpline,0);
            T val2 = (*it).second;
            cout << (*it).first << ": " << val1 << " " << val2 << " " << fabs(val1-val2) << " "
                 << fabs(val1-val2)/fabs(val1) << endl;
        }
    }


    return 0;
}

