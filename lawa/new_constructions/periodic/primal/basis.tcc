#ifndef LAWA_CONSTRUCTIONS_PERIODIC_PRIMAL_BASIS_TCC
#define LAWA_CONSTRUCTIONS_PERIODIC_PRIMAL_BASIS_TCC 1

#include <cassert>

namespace lawa {

template <typename T>
Basis<T,Primal,Periodic,CDF>::Basis(int _d, int _d_, int j)
    : d(_d), d_(_d_), j0(j), mra(d,d_,j), mra_(d,d_,j), 
      psi(d,d_), M1(psi), _j(j)
{
}

template <typename T>
int
Basis<T,Primal,Periodic,CDF>::level() const
{
    return _j;
}

template <typename T>
void
Basis<T,Primal,Periodic,CDF>::setLevel(int j) const
{
    assert(j>=j0);
    _j = j;
}

template <typename T>
const BasisFunction<T,Primal,Periodic,CDF> &
Basis<T,Primal,Periodic,CDF>::generator(XType xtype) const
{
    if (xtype==XBSpline) {
        return mra.phi;
    } else {
        return psi;
    }
}

template <typename T>
int
Basis<T,Primal,Periodic,CDF>::cardJ(int j) const
{
    assert(j>=j0);

    return pow2i<T>(j);
}

template <typename T>
Range<Integer>
Basis<T,Primal,Periodic,CDF>::rangeJ(int j) const
{
    assert(j>=j0);

    return Range<Integer>(1,pow2i<T>(j));
}

} // namespace lawa

#endif // LAWA_CONSTRUCTIONS_PERIODIC_PRIMAL_BASIS_TCC
