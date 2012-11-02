/* TEST LOCAL RECONSTRUCTION
 *
 *  This examples calculates the local single reconstruction for wavelet interval bases. We consider
 *  Dijkema wavelet as well as L2-orthonormal multiwavelets.
 *
 */

#include <iostream>
#include <fstream>
#include <lawa/lawa.h>
#include <lawa/methods/adaptive/datastructures/treecoefficients1d.h>

using namespace std;
using namespace lawa;

/// Several typedefs for notational convenience.

///  Typedefs for Flens data types:
typedef double T;
typedef flens::DenseVector<flens::Array<T> >                        DenseVectorT;
typedef flens::GeMatrix<flens::FullStorage<T, cxxblas::ColMajor> >  DenseMatrixT;

///  Typedefs for problem components:
///  Wavelet basis over an interval
typedef Basis<T, Orthogonal, R, Multi>                       PrimalBasis;
typedef PrimalBasis::RefinementBasis                                RefinementBasis;

///  Wavelet integrals
typedef Integral<Gauss,RefinementBasis,RefinementBasis>             IntegralRefinentBasis;
typedef IntegralF<Gauss,PrimalBasis>                                IntegralFBasis;

typedef CoefficientsByLevel<T>::const_it                            const_coeffbylevel_it;
typedef IndexSet<Index1D>::const_iterator                           const_set1d_it;
typedef Coefficients<Lexicographical,T,Index1D>::const_iterator     const_coeff1d_it;

int K = 10;

T
f(T x) { return std::exp(-10.L*(x-0.5L)*(x-0.5)); }

void
constructRandomTree(const PrimalBasis &basis, int J, TreeCoefficients1D<T> &LambdaTree);

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
    int j  = j0+5;
    int J  = atoi(argv[3]);

    /// Basis initialization, using Dirichlet boundary conditions
    PrimalBasis basis(d, j0);           // For L2_orthonormal and special MW bases
    RefinementBasis &refinementbasis = basis.refinementbasis;

    /// Integral initialization
    DenseVectorT    singPts;
    Function<T>     f_fct(f,singPts);
    IntegralFBasis  integralf(f_fct,basis);
    integralf.quadrature.setOrder(40);

    /// Local refinement initialization
    LocalRefinement<PrimalBasis> LocalRefine(basis);

    /// Computing a vector of (multi-)scaling coefficients and a vector of (multi-)wavelet coefficients
    CoefficientsByLevel<T> u_scaling1, u_wavelet1;
    for (int k=-2*K; k<=2*K; ++k) {
        u_scaling1.map.operator[](k) = integralf(j,k,XBSpline,0);
    }
    for (int k=-K; k<=K; ++k) {
        u_wavelet1.map.operator[](k) = integralf(j,k,XWavelet,0);
    }

    /// In case of multiscaling functions, we have to represent them first in a B-spline basis
    /// (refinement basis).
    CoefficientsByLevel<T> u_bspline1;
    int refinement_j_bspline = 0;
    if (PrimalBasis::Cons==Multi && d>1) {
        cout << "Reconstruct multi scaling function..." << endl;
        LocalRefine.reconstructOnlyMultiScaling(u_scaling1, j, u_bspline1, refinement_j_bspline);
        std::cout << "refinement_j_bspline = " << refinement_j_bspline << std::endl;
    }
    else {
        refinement_j_bspline = j;
        u_bspline1 = u_scaling1;
    }

    /// Now we compute the common representation of the b-spline coefficient vector and the (multi)-
    /// wavelet coefficient vector.
    CoefficientsByLevel<T> u_loc_single1;
    int refinement_j = 0;
    Timer time;
    time.start();
    LocalRefine.reconstruct(u_bspline1, refinement_j_bspline, u_wavelet1, j, u_loc_single1, refinement_j);
    time.stop();
    cout << "Local reconstruction took: " << time.elapsed() << endl;

    ofstream plotfile("test.txt");
    plotfile.precision(16);
    for (T x=-10; x<=10; x+=pow2i<T>(-8-j)) {
        T val1 = 0.L, val2 = 0.L, val3 = 0.L;
        // no refinement
        for (const_coeffbylevel_it it=u_scaling1.map.begin(); it!=u_scaling1.map.end(); ++it) {
            val1 += (*it).second * basis.generator(XBSpline).operator()(x,j,(*it).first,0);
        }
        for (const_coeffbylevel_it it=u_wavelet1.map.begin(); it!=u_wavelet1.map.end(); ++it) {
            val1 += (*it).second * basis.generator(XWavelet).operator()(x,j,(*it).first,0);
        }
        // from multiscaling to bspline representation
        for (const_coeffbylevel_it it=u_bspline1.map.begin(); it!=u_bspline1.map.end(); ++it) {
            val2 += (*it).second * refinementbasis.generator(XBSpline).operator()(x,refinement_j_bspline,(*it).first,0);
        }
        for (const_coeffbylevel_it it=u_wavelet1.map.begin(); it!=u_wavelet1.map.end(); ++it) {
            val2 += (*it).second * basis.generator(XWavelet).operator()(x,j,(*it).first,0);
        }
        // single scale representation
        for (const_coeffbylevel_it it=u_loc_single1.map.begin(); it!=u_loc_single1.map.end(); ++it) {
            val3 += (*it).second * refinementbasis.generator(XBSpline).operator()(x,refinement_j,(*it).first,0);
        }
        plotfile << x << " " << f(x) << " " << val1 << " " << val2 << " " << val3 << endl;
    }


    /// We construct a random tree on which we test the transformation to a local scaling function
    /// representation
    /*
    TreeCoefficients1D<T> u_tree(4096,basis.j0);
    Coefficients<Lexicographical,T,Index1D> u, u_loc_single;
    constructRandomTree(basis, J, u_tree);
    fromTreeCoefficientsToCoefficients(u_tree,u);
    IndexSet<Index1D> supp_u;
    supp_u = supp(u);
    for (const_set1d_it it=supp_u.begin(); it!=supp_u.end(); ++it) {
        u[*it] = integralf((*it).j,(*it).k,(*it).xtype,0);
    }

    /// The vector u contains the multilevel representation of a coefficient vector. We transform
    /// it to the local single scale representation.
    time.start();
    LocalRefine.reconstruct(u, j0, u_loc_single);
    time.stop();
    cout << "Local reconstruction took " << time.elapsed() << endl;


    /// We plot the result. Note that the function f is only approximated if L2-orthonormal multi
    /// wavelets are used!!
    ofstream plotfile2("test2.txt");
    plotfile2.precision(16);
    T max_error1=0.L, max_error2=0.L;
    for (T x=0.; x<1.; x+=pow2i<T>(-6-J)) {
        T val1=0.L, val2=0.L;
        for (const_coeff1d_it it=u.begin(); it!=u.end(); ++it) {
            val1 += (*it).second *
                    basis.generator((*it).first.xtype).operator()(x,(*it).first.j,(*it).first.k,0);
        }
        for (const_coeff1d_it it=u_loc_single.begin(); it!=u_loc_single.end(); ++it) {
            val2 += (*it).second *
                    refinementbasis.mra.phi.operator()(x,(*it).first.j,(*it).first.k,0);
        }
        max_error1 = std::max(max_error1,fabs(f(x)-val1));
        max_error2 = std::max(max_error2,fabs(val1-val2));
        plotfile2 << x << " " << f(x) << " " << val1 << " " << val2 << endl;
    }
    plotfile2.close();
    cout << u.size() << " " << max_error1 << " " << max_error2 << endl;

    /// Finally, we visualize the corresponding index sets.
    plotCoeff(u, basis, "coeff_multi_scale", false, true);
    plotCoeff(u_loc_single, basis, "coeff_local_single_scale", true, true);
    */
    return 0;
}

void
constructRandomTree(const PrimalBasis &basis, int J, TreeCoefficients1D<T> &LambdaTree)
{
    /*
    for (int k=basis.mra.rangeI(basis.j0).firstIndex(); k<=basis.mra.rangeI(basis.j0).lastIndex(); ++k) {
        LambdaTree[0].map.operator[](k) = 0.;
    }
    for (int j=basis.j0; j<=J; ++j) {
        for (int k=basis.rangeJ(j).firstIndex(); k<=basis.rangeJ(j).lastIndex(); ++k) {
            LambdaTree.bylevel[j-basis.j0+1].map.operator[](k) = 0.;
        }
    }
    */

    for (int k=-2*K; k<=2*K; ++k) {
        LambdaTree[0].map.operator[](k) = 0.;
    }
    for (int j=basis.j0; j<=J; ++j) {
        int random_k1 = rand() % (2*K+1) + 1;
        LambdaTree.bylevel[j-basis.j0+1].map.operator[](random_k1) = 0.;
        int random_k2 = rand() % (2*K+1) + 1;
        LambdaTree.bylevel[j-basis.j0+1].map.operator[](random_k2) = 0.;
    }
    for (int i=J-basis.j0+1; i>=2; --i) {
        CoefficientsByLevel<T> *currentlevel;
        currentlevel = &(LambdaTree.bylevel[i]);
        int j=basis.j0+i-1;
        for (const_coeffbylevel_it it=(*currentlevel).map.begin(); it!=(*currentlevel).map.end(); ++it) {
            long k = (*it).first;
            long k_first = (int)k / 2 - 30;
            long k_last  = (int)k / 2 + 30;
            for (int k1=k_first; k1<=k_last; ++k1) {
                if (overlap(basis.psi.support(j,k),basis.psi.support(j-1,k1))>0) {
                    LambdaTree.bylevel[i-1].map.operator[](k1) = 0.;
                }
            }
        }
    }
}
