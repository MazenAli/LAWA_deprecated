#include <cmath>

namespace lawa {

template <typename T, typename Basis2D>
LeftNormPreconditioner2D<T,Basis2D>::LeftNormPreconditioner2D(const Basis2D &basis, T s)
    : _basis(basis), _s(s),
      _integral(basis.second, basis.second)
{
}

template <typename T, typename Basis2D>
T
LeftNormPreconditioner2D<T,Basis2D>::operator()(XType xtype1, int j1, int k1,
                                                XType xtype2, int j2, int k2) const
{
    T value = _integral(xtype2,j2,k2,0, xtype2,j2,k2,0);

    if (_s==2.) {
        return 1./std::sqrt(value + _integral(xtype2,j2,k2,1, xtype2,j2,k2,1));
    } else {
        return 1./std::sqrt(value + std::pow(2.,_s*j2));
    }
}

template <typename T, typename Basis2D>
T
LeftNormPreconditioner2D<T,Basis2D>::operator()(const Index2D &index) const
{
    return this->operator()(index.index1.xtype, index.index1.j, index.index1.k,
                            index.index2.xtype, index.index2.j, index.index2.k);
}

}   // namespace lawa