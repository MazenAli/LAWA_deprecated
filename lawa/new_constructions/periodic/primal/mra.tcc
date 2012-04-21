#ifndef LAWA_CONSTRUCTIONS_PERIODIC_PRIMAL_MRA_TCC
#define LAWA_CONSTRUCTIONS_PERIODIC_PRIMAL_MRA_TCC 1

#include <cassert>

namespace lawa {

template <typename T>
MRA<T,Primal,Periodic,CDF>::MRA(int _d, int _d_, int j)
    : d(_d), d_(_d_), j0(j), phi(d), M0(phi), _j(j)
{
}

template <typename T>
int
MRA<T,Primal,Periodic,CDF>::level() const
{
    return _j;
}

template <typename T>
void
MRA<T,Primal,Periodic,CDF>::setLevel(int j) const
{
    assert(j>=0);
    _j = j;
}

template <typename T>
int
MRA<T,Primal,Periodic,CDF>::cardI(int j) const
{
    assert(j>=j0);   
    return pow2i<T>(j);
}

template <typename T>
Range<Integer>
MRA<T,Primal,Periodic,CDF>::rangeI(int j) const
{
    assert(j>=j0);
    return Range<Integer>(1,pow2i<T>(j));
}

} // namespace lawa

#endif // LAWA_CONSTRUCTIONS_PERIODIC_PRIMAL_MRA_TCC
