#include <iostream>
#include <fstream>
#include <lawa/lawa.h>
#include <applications/finance/operators/cgmyoperator1d.h>
#include <applications/finance/righthandsides/righthandsides.h>

using namespace std;
using namespace lawa;

typedef long double T;

typedef Basis<T, Orthogonal, Interval, MultiRefinement>                       PrimalBasis;

struct CGMYKernel {

    CGMYKernel(const Kernel<T,CGMY> &_kernel) : kernel(_kernel) { };

    T
    operator()(T x) const { return kernel.SecondTailIntegral(-x); }

    const Kernel<T,CGMY> &kernel;
};

int main()
{
    cout.precision(16);
    /// wavelet basis parameters:
    int d = 2;          // mandatory for this test case!!!
    int j0 = 0;         // minimal level

    int j1 = 3, k1 = 4; XType xtype1 = XBSpline;
    int j2 = 3, k2 = 3; XType xtype2 = XBSpline;

    ProcessParameters1D<T,CGMY> processparameters(0., 1., 2.4, 4.8, 1.1);
    Kernel<T,CGMY> kernel(processparameters);
    CGMYKernel cgmykernel(kernel);

    PrimalBasis basis(d, j0);
    basis.enforceBoundaryCondition<DirichletBC>();

    GeMatrix<FullStorage<T,ColMajor> > varphi_row_deltas, varphi_col_deltas;
    varphi_row_deltas = computeDeltas<T,PrimalBasis>(basis,j1,k1,xtype1);
    varphi_col_deltas = computeDeltas<T,PrimalBasis>(basis,j2,k2,xtype2);

    T part1=(T)0., part2=(T)0.;
    for (int lambda=varphi_row_deltas.rows().firstIndex(); lambda<=varphi_row_deltas.rows().lastIndex(); ++lambda) {
        T x = varphi_row_deltas(lambda,1);

        part1 += varphi_row_deltas(lambda,2)*basis.generator(xtype2)(x,j2,k2,0)*kernel.c3;

        for (int mu=varphi_col_deltas.rows().firstIndex(); mu<=varphi_col_deltas.rows().lastIndex(); ++mu) {
            T y = varphi_col_deltas(mu,1);
            T c = varphi_col_deltas(mu,2)*varphi_row_deltas(lambda,2);

            if (fabs(x-y)>1e-10)  {
                T value_tailintegral=kernel.ForthTailIntegral(y-x);
                if (y-x>0)  part2 += c * (value_tailintegral - kernel.constants[2]);
                else        part2 += c * (value_tailintegral - kernel.constants[3]);
            }
        }
    }
    T integral_value1 = -(part1 + part2);

    SingularIntegral<CGMYKernel,PrimalBasis,PrimalBasis> singularIntegral(cgmykernel,basis,basis);
    int order = 15, n = 10;
    T sigma = 0.1, mu = 0.3, omega = 0.01;

    for (n=15; n<=40; ++n) {
        singularIntegral.singularquadrature.setParameters(order, n, sigma, mu, omega);
        T integral_value2 = singularIntegral(j1,k1,xtype1,1, j2,k2,xtype2,1);

        cout << "integral_value1 = " << integral_value1 << ", integral_value2 = " << integral_value2
             << ", diff = " << integral_value1 - integral_value2 << endl;

    }
    return 0;

}
