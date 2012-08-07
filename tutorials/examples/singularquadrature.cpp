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
///     Primal Basis over an interval, using Dijkema construction
typedef Basis<T, Orthogonal, Interval, Multi>                         PrimalBasis;

struct LogKernel {

    LogKernel(void) { };

    T
    operator()(T x) const { return log(fabs(x)); }
};

int main()
{
    /// wavelet basis parameters:
    int d = 1;          // (d,d_)-wavelets
    int d_ = 1;
    int j0 = 0;         // minimal level
    int J = 5;          // maximal level

    /// Basis initialization, using Dirichlet boundary conditions
    //PrimalBasis basis(d, d_, j0);
    PrimalBasis basis(d, j0);

    if (d>=2) basis.enforceBoundaryCondition<DirichletBC>();

    LogKernel logkernel;

    SingularIntegral<LogKernel,PrimalBasis,PrimalBasis> singularIntegral(logkernel,basis,basis);

    int order = 2;
    int n = 10;
    T sigma = 0.2;
    T mu = 0.3;
    T omega = 0.01;


    cout.precision(16);
    singularIntegral.singularquadrature.setParameters(order, n, sigma, mu, omega);
    cout << "order = " << order << ", n = " << n << ", sigma = " << sigma << ", mu = " << mu << ": "
         << singularIntegral(0, 0, XBSpline, 0, 0, 0, XBSpline, 0)+1.5 << endl;

    mu = 0.5;
    singularIntegral.singularquadrature.setParameters(order, n, sigma, mu, omega);
    cout << "order = " << order << ", n = " << n << ", sigma = " << sigma << ", mu = " << mu << ": "
         << singularIntegral(0, 0, XBSpline, 0, 0, 0, XBSpline, 0)+1.5 << endl;

    mu = 0.5;
    order = 1;
    singularIntegral.singularquadrature.setParameters(order, n, sigma, mu, omega);
    cout << "order = " << order << ", n = " << n << ", sigma = " << sigma << ", mu = " << mu << ": "
         << singularIntegral(0, 0, XBSpline, 0, 0, 0, XBSpline, 0)+1.5 << endl;

    sigma = 0.1;
    singularIntegral.singularquadrature.setParameters(order, n, sigma, mu, omega);
    cout << "order = " << order << ", n = " << n << ", sigma = " << sigma << ", mu = " << mu << ": "
         << singularIntegral(0, 0, XBSpline, 0, 0, 0, XBSpline, 0)+1.5 << endl;

    sigma = 0.3;
    singularIntegral.singularquadrature.setParameters(order, n, sigma, mu, omega);
    cout << "order = " << order << ", n = " << n << ", sigma = " << sigma << ", mu = " << mu << ": "
         << singularIntegral(0, 0, XBSpline, 0, 0, 0, XBSpline, 0)+1.5 << endl;

    sigma = 0.2;
    n = 15;
    singularIntegral.singularquadrature.setParameters(order, n, sigma, mu, omega);
    cout << "order = " << order << ", n = " << n << ", sigma = " << sigma << ", mu = " << mu << ": "
         << singularIntegral(0, 0, XBSpline, 0, 0, 0, XBSpline, 0)+1.5 << endl;

    n = 20;
    singularIntegral.singularquadrature.setParameters(order, n, sigma, mu, omega);
    cout << "order = " << order << ", n = " << n << ", sigma = " << sigma << ", mu = " << mu << ": "
         << singularIntegral(0, 0, XBSpline, 0, 0, 0, XBSpline, 0)+1.5 << endl;

    return 0;
}
