#include <cmath>

namespace lawa {

template <typename T, typename Basis>
H1NormPreconditioner1D<T,Basis>::H1NormPreconditioner1D(const Basis &basis)
    : _integral(basis,basis)
{
}

template <typename T, typename Basis>
T
H1NormPreconditioner1D<T,Basis>::operator()(XType xtype, int j, int k) const
{
    return 1./std::sqrt(_integral(xtype,j,k,0, xtype,j,k,0)
                      + _integral(xtype,j,k,1, xtype,j,k,1));
}

template <typename T, typename Basis>
T
H1NormPreconditioner1D<T,Basis>::operator()(const Index1D &index) const
{
    return this->operator()(index.xtype,index.j,index.k);
}

}   // namespace lawa
