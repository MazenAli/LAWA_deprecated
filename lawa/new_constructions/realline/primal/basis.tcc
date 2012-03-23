#ifndef LAWA_CONSTRUCTIONS_REALLINE_PRIMAL_BASIS_TCC
#define LAWA_CONSTRUCTIONS_REALLINE_PRIMAL_BASIS_TCC 1

namespace lawa {

template <typename T>
Basis<T,Primal,R,CDF>::Basis(int _d, int _d_, int j)
    : d(_d), d_(_d_), j0(j), mra(d), mra_(d,d_), psi(d,d_), M1(psi)
{
}

template <typename T>
int
Basis<T,Primal,R,CDF>::level() const
{
    return _j;
}

template <typename T>
void
Basis<T,Primal,R,CDF>::setLevel(int j) const
{
    _j = j;
}

template <typename T>
const BasisFunction<T,Primal,R,CDF> &
Basis<T,Primal,R,CDF>::generator(XType xtype) const
{
    if (xtype==XBSpline) {
        return mra.phi;
    } else {
        return psi;
    }
}

} // namespace lawa

#endif // LAWA_CONSTRUCTIONS_REALLINE_PRIMAL_BASIS_TCC
