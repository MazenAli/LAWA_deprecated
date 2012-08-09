#include <iostream>
#include <fstream>
#include <lawa/lawa.h>

using namespace std;
using namespace lawa;

/// Several typedefs for notational convenience.

///  Typedefs for Flens data types:
typedef double T;
typedef flens::GeMatrix<flens::FullStorage<T, cxxblas::ColMajor> >  FullColMatrixT;
typedef flens::SparseGeMatrix<flens::CRS<T,flens::CRS_General> >    SparseMatrixT;
typedef flens::DiagonalMatrix<T>                                    DiagonalMatrixT;
typedef flens::DenseVector<flens::Array<T> >                        DenseVectorT;

///  Typedefs for problem components:
typedef Basis<T, Orthogonal, Interval, Multi>                         PrimalBasis;

struct LogKernel {

    LogKernel(void) { };

    T
    operator()(T x) const { return log(fabs(x)); }
};

int main()
{
    cout.precision(16);
    /// wavelet basis parameters:
    int j0 = 0;         // minimal level

    LogKernel logkernel;
    Poly<T> p1(2), p2(2);
    SingularIntegralPP<LogKernel,Poly<T>,Poly<T> > singularIntegralPP(logkernel,p1,p2);

    int order = 7, n = 20;
    T sigma = 0.2, mu = 0.5, omega = 0.01;
    T ref_val = 0., approx = 0.;
    singularIntegralPP.singularquadrature.setParameters(order,n,sigma,mu,omega);

    // Non-singular integral, non-quadratic domain
    ref_val = (-340. + 10824.*log(2.) - 1296.*log(3.) + 1225*log(5.) - 3969*log(7.))/32768.;
    approx  = singularIntegralPP(0.,0.125,0.75,1.);
    cout << "approx = " << approx << ", ref_val = " << ref_val << ", diff = " << approx-ref_val << endl;


    // Corner singularity in (a1,b2), quadratic domain
    ref_val = (-217. + log(256.))/4096.;
    approx  = singularIntegralPP(0.75,1.,0.5,0.75);
    cout << "approx = " << approx << ", ref_val = " << ref_val << ", diff = " << approx-ref_val << endl;

    // Corner singularity in (a1,b2), non-quadratic domain 1
    ref_val = (-3380-6664*log(2.) - 1225*log(5.) + 3969*log(7.))/32768.;
    approx  = singularIntegralPP(0.75,1.,0.125,0.75);
    cout << "approx = " << approx << ", ref_val = " << ref_val << ", diff = " << approx-ref_val << endl;

    // Corner singularity in (a1,b2), non-quadratic domain 2
    ref_val = (-1225. + 23820.*log(2.) - 7938.*log(7.))/65536.;
    approx  = singularIntegralPP(0.125,1.,0.,0.125);
    cout << "approx = " << approx << ", ref_val = " << ref_val << ", diff = " << approx-ref_val << endl;


    // Corner singularity in (b1,a2), quadratic domain
    ref_val = (-97.-52.*log(2.))/65536.;
    approx  = singularIntegralPP(0.125,0.25,0.25,0.375);
    cout << "approx = " << approx << ", ref_val = " << ref_val << ", diff = " << approx-ref_val << endl;

    // Corner singularity in (b1,a2) non-quadratic domain
    ref_val = -3.*(272. + 1560.*log(2.) + 1200.*log(3.) - 1323*log(7.))/32768.;
    approx  = singularIntegralPP(0.125,0.25,0.25,1.);
    cout << "approx = " << approx << ", ref_val = " << ref_val << ", diff = " << approx-ref_val << endl;

    // Corner singularity in (b1,a2) non-quadratic domain
    ref_val = (-3065. + 5040.*log(8./5.) - 15972.*log(2.) + 4608*log(3.) + 2590.*log(5.))/65536.;
    approx  = singularIntegralPP(0.125,0.75,0.75,0.875);
    cout << "approx = " << approx << ", ref_val = " << ref_val << ", diff = " << approx-ref_val << endl;


    // Diagonal singularity in non-quadratic domain
    ref_val = (-10225. - 10196*log(2.) + 4608.*log(3.))/65536.;
    approx  = singularIntegralPP(0.125,0.75,0.25,0.875);
    cout << "approx = " << approx << ", ref_val = " << ref_val << ", diff = " << approx-ref_val << endl;

    return 0;
}


/*
 *
 *
    /// Basis initialization,
    PrimalBasis basis1(1, j0);
    PrimalBasis basis2(2, j0);

    /// Using Dirichlet boundary conditions for higher order bases
    basis2.enforceBoundaryCondition<DirichletBC>();



    SingularIntegral<LogKernel,PrimalBasis,PrimalBasis> singularIntegral1(logkernel,basis1,basis1);
    SingularIntegral<LogKernel,PrimalBasis,PrimalBasis> singularIntegral2(logkernel,basis2,basis2);

    /// Test singular quadrature for constant functions
    cout << "*** Testing constant basis functions ***" << endl;
    int order1 = 1, n1 = 20;
    T sigma1 = 0.2, mu1 = 0.5, omega1 = 0.01;
    T ref_val1 = 0., approx1 = 0.;

    singularIntegral1.singularquadrature.setParameters(order1, n1, sigma1, mu1, omega1);

    ref_val1 = -1.5;
    approx1  = singularIntegral1(0, 0, XBSpline, 0, 0, 0, XBSpline, 0);
    cout << "order = " << order1 << ", n = " << n1 << ", sigma = " << sigma1 << ", mu = " << mu1 << ": "
         << approx1 << " " << ref_val1 << " " << approx1 - ref_val1 << endl;

    ref_val1 = 2.*(1./8.)*(-3-2*log(2.));
    approx1  = singularIntegral1(1, 0, XBSpline, 0, 1, 0, XBSpline, 0);
    cout << "order = " << order1 << ", n = " << n1 << ", sigma = " << sigma1 << ", mu = " << mu1 << ": "
         << approx1 << " " << ref_val1 << " " << approx1 - ref_val1 << endl;

    ref_val1 = 2.*(1./8.)*(-3+log(4.));
    approx1  = singularIntegral1(1, 0, XBSpline, 0, 1, 1, XBSpline, 0);
    cout << "order = " << order1 << ", n = " << n1 << ", sigma = " << sigma1 << ", mu = " << mu1 << ": "
         << approx1 << " " << ref_val1 << " " << approx1 - ref_val1 << endl;

    ref_val1 = 2.*(1./8.)*(-3+log(4.));
    approx1  = singularIntegral1(1, 1, XBSpline, 0, 1, 0, XBSpline, 0);
    cout << "order = " << order1 << ", n = " << n1 << ", sigma = " << sigma1 << ", mu = " << mu1 << ": "
          << approx1 << " " << ref_val1 << " " << approx1 - ref_val1 << endl;

    cout << endl << endl;

    /// Test singular quadrature for linear functions
    cout << "*** Testing linear basis functions ***" << endl;
    int order2 = 2, n2 = 20;
    T sigma2 = 0.2, mu2 = 0.5, omega2 = 0.01;
    T ref_val2 = 0., approx2 = 0.;

    //ref_val2 = (1./64.)*(-100.+28.*log(2.)-3.*log(16.));
    ref_val2 = (-45. + 616.*log(2.) - 756.*log(3.) + 320*log(4.))/1024.;
    //for (n2=5; n2<=40; n2+=5) {

    for (order2=5; order2<=120; order2+=5) {
        singularIntegral2.singularquadrature.setParameters(order2, n2, sigma2, mu2, omega2);
        approx2  = singularIntegral2(0, 0, XBSpline, 0, 0, 0, XBSpline, 0);
        cout << "order = " << order2 << ", n = " << n2 << ", sigma = " << sigma2 << ", mu = " << mu2 << ": "
             << approx2 << " " << ref_val2 << " " << approx2 - ref_val2 << endl;
    }
 */
