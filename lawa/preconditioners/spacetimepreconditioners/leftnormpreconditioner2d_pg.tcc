#include <cmath>

namespace lawa {

template <typename T, typename TrialBasis, typename TestBasis>
LeftNormPreconditioner2D_PG<T,TrialBasis, TestBasis>::LeftNormPreconditioner2D_PG(const TrialBasis &trialbasis, const TestBasis &testbasis, T s)
    : _trialbasis(trialbasis), _testbasis(testbasis), _s(s),
      _integral(testbasis.second, trialbasis.second)
{
}

template <typename T, typename TrialBasis, typename TestBasis>
T
LeftNormPreconditioner2D_PG<T,TrialBasis, TestBasis>::operator()(XType /*xtype1*/, int /*j1*/, int /*k1*/,
                                                XType xtype2, int j2, int k2) const
{
    T value = _integral(j2,k2,xtype2,0,j2,k2,xtype2,0);

    if (_s==2.) {
        // Calculate H1-Norm of Basis Function using Integrals
        return 1./std::sqrt(value + _integral(j2,k2,xtype2,1,j2,k2,xtype2,1));
    } else {
        // Calculate H1-Norm of Basis Function using Scaling of the L2-norm (assumed to be equivalent to 1)
        return 1./std::sqrt(value + std::pow(2.,_s*j2));
    }
}

template <typename T, typename TrialBasis, typename TestBasis>
T
LeftNormPreconditioner2D_PG<T,TrialBasis, TestBasis>::operator()(const Index2D &index) const
{
    return this->operator()(index.index1.xtype, index.index1.j, index.index1.k,
                            index.index2.xtype, index.index2.j, index.index2.k);
}

}   // namespace lawa

