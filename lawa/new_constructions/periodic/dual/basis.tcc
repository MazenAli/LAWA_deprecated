#ifndef LAWA_CONSTRUCTIONS_PERIODIC_DUAL_BASIS_TCC
#define LAWA_CONSTRUCTIONS_PERIODIC_DUAL_BASIS_TCC 1

#include <cassert>

namespace lawa {

template <typename T>
Basis<T,Dual,Periodic,CDF>::Basis(unsigned int _d, unsigned int _d_, int j)
    : d(_d), d_(_d_), j0(j), mra(d,d_,j), mra_(d,d_,j), 
      psi_(d,d_), M1_(psi_), _j(j)
{
}

template <typename T>
int
Basis<T,Dual,Periodic,CDF>::level() const
{
    return _j;
}

template <typename T>
void
Basis<T,Dual,Periodic,CDF>::setLevel(int j) const
{
    assert(j>=j0);
    _j = j;
}

template <typename T>
const BasisFunction<T,Dual,Periodic,CDF> &
Basis<T,Dual,Periodic,CDF>::generator(XType xtype) const
{
    if (xtype==XBSpline) {
        return mra.phi_;
    } else {
        return psi_;
    }
}

template <typename T>
int
Basis<T,Dual,Periodic,CDF>::cardJ_(int j) const
{
    assert(j>=j0);

    return pow2i<T>(j);
}

template <typename T>
Range<Integer>
Basis<T,Dual,Periodic,CDF>::rangeJ_(int j) const
{
    assert(j>=j0);

    return Range<Integer>(0,pow2i<T>(j)-1);
}

} // namespace lawa

#endif // LAWA_CONSTRUCTIONS_PERIODIC_DUAL_BASIS_TCC
