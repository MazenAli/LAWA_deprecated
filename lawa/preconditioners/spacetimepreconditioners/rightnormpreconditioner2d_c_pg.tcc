#include <cmath>

namespace lawa {

template <typename T, typename TrialBasis, typename TestBasis>
RightNormPreconditioner2D_c_PG<T,TrialBasis, TestBasis>::RightNormPreconditioner2D_c_PG(const TrialBasis& trialbasis, const TestBasis &testbasis, T s)
    : _s(s), _integral_t(testbasis.first,trialbasis.first), 
             _integral_x(testbasis.second,trialbasis.second)
{
}

template <typename T, typename TrialBasis, typename TestBasis>
T
RightNormPreconditioner2D_c_PG<T,TrialBasis, TestBasis>::operator()(XType xtype1, int j1, int k1,
                                                 XType xtype2, int j2, int k2) const
{
    T value_t    = _integral_t(j1,k1,xtype1,0, j1,k1,xtype1,0);
    //T dd_value_t = _integral_t(j1,k1,xtype1,1, j1,k1,xtype1,1);
    
    T value_x    = _integral_x(j2,k2,xtype2,0, j2,k2,xtype2,0);
    T dd_value_x = _integral_x(j2,k2,xtype2,1, j2,k2,xtype2,1);

    if (_s==2.) {
        return 1./std::sqrt( value_t*(value_x+dd_value_x) + pow2i<T>(2*j1)*value_t/(value_x + dd_value_x));
    }
    else {
        return 1./std::sqrt((value_x+std::pow(2.,_s*j2)) + std::pow(2.,-_s*j1)*value_t*std::pow(2.,-_s*j2));
    }
}

template <typename T, typename TrialBasis, typename TestBasis>
T
RightNormPreconditioner2D_c_PG<T,TrialBasis, TestBasis>::operator()(const Index2D &index) const
{
    return this->operator()(index.index1.xtype, index.index1.j, index.index1.k,
                            index.index2.xtype, index.index2.j, index.index2.k);
}

}   // namespace lawa
