#ifndef LAWA_CONSTRUCTIONS_PERIODIC_DUAL_MRA_TCC
#define LAWA_CONSTRUCTIONS_PERIODIC_DUAL_MRA_TCC 1

#include <cassert>

namespace lawa {

template <typename T>
MRA<T,Dual,Periodic,CDF>::MRA(int _d, int _d_, int j)
    : d(_d), d_(_d_), j0(j),
      phi_(d,d_), M0_(phi_), _j(j)
{
}

template <typename T>
int
MRA<T,Dual,Periodic,CDF>::level() const
{
    return _j;
}

template <typename T>
void
MRA<T,Dual,Periodic,CDF>::setLevel(int j) const
{
    assert(j>=0);
    _j = j;
}

template <typename T>
int
MRA<T,Dual,Periodic,CDF>::cardI_(int j) const
{
    assert(j>=j0);   
    return pow2i<T>(j);
}

template <typename T>
Range<Integer>
MRA<T,Dual,Periodic,CDF>::rangeI_(int j) const
{
    assert(j>=j0);
    assert(0);
    return Range<Integer>(0,pow2i<T>(j)-1);
}

} // namespace lawa

#endif // LAWA_CONSTRUCTIONS_PERIODIC_DUAL_MRA_TCC
