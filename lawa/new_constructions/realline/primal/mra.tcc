#ifndef LAWA_CONSTRUCTIONS_REALLINE_PRIMAL_MRA_TCC
#define LAWA_CONSTRUCTIONS_REALLINE_PRIMAL_MRA_TCC 1

namespace lawa {

template <typename T>
MRA<T,Primal,R,CDF>::MRA(int _d, int j)
    : d(_d), j0(j), phi(d), M0(phi)
{
}

template <typename T>
int
MRA<T,Primal,R,CDF>::level() const
{
    return _j;
}

template <typename T>
void
MRA<T,Primal,R,CDF>::setLevel(int j) const
{
    _j = j;
}

} // namespace lawa

#endif // LAWA_CONSTRUCTIONS_REALLINE_PRIMAL_MRA_TCC
