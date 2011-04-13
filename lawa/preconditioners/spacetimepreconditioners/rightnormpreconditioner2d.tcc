#include <cmath>

namespace lawa {

template <typename T, typename Basis2D>
RightNormPreconditioner2D<T,Basis2D>::RightNormPreconditioner2D(const Basis2D &basis, T s)
    : _s(s), _integral_t(basis.first,basis.first), 
             _integral_x(basis.second,basis.second)
{
}

template <typename T, typename Basis2D>
T
RightNormPreconditioner2D<T,Basis2D>::operator()(XType xtype1, int j1, int k1,
                                                 XType xtype2, int j2, int k2) const
{
    T value_t    = _integral_t(xtype1,j1,k1,0, xtype1,j1,k1,0);
    T dd_value_t = _integral_t(xtype1,j1,k1,1, xtype1,j1,k1,1);
    
    T value_x    = _integral_x(xtype2,j2,k2,0, xtype2,j2,k2,0);
    T dd_value_x = _integral_x(xtype2,j2,k2,1, xtype2,j2,k2,1);

    if (_s==2.) {
        return 1./std::sqrt( (value_x+dd_value_x) + (value_t+dd_value_t)*pow2i<T>(-2*j2));
    }
    else {
        return 1./std::sqrt((value_x+std::pow(2.,_s*j2)) + (value_t+dd_value_t)*std::pow(2.,-_s*j2));
    }
}

template <typename T, typename Basis2D>
T
RightNormPreconditioner2D<T,Basis2D>::operator()(const Index2D &index) const
{
    return this->operator()(index.index1.xtype, index.index1.j, index.index1.k,
                            index.index2.xtype, index.index2.j, index.index2.k);
}

}   // namespace lawa