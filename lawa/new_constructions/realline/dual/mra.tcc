#ifndef LAWA_CONSTRUCTIONS_REALLINE_DUAL_MRA_TCC
#define LAWA_CONSTRUCTIONS_REALLINE_DUAL_MRA_TCC 1

namespace lawa {

template <typename T>
MRA<T,Dual,R,CDF>::MRA(int _d, int _d_, int j)
    : d(_d), d_(_d_), j0(j), phi_(d,d_), M0_(phi_)
{
}

template <typename T>
int
MRA<T,Dual,R,CDF>::level() const
{
    return _j;
}

template <typename T>
void
MRA<T,Dual,R,CDF>::setLevel(int j) const
{
    _j = j;
}

} // namespace lawa

#endif // LAWA_CONSTRUCTIONS_REALLINE_DUAL_MRA_TCC
