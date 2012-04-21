#ifndef LAWA_CONSTRUCTIONS_REALLINE_DUAL_BASIS_TCC
#define LAWA_CONSTRUCTIONS_REALLINE_DUAL_BASIS_TCC 1

namespace lawa {

template <typename T>
Basis<T,Dual,R,CDF>::Basis(int _d, int _d_, int j)
    : d(_d), d_(_d_), j0(j), mra(d), mra_(d,d_), psi_(d,d_), M1_(psi_)
{
}

template <typename T>
int
Basis<T,Dual,R,CDF>::level() const
{
    return _j;
}

template <typename T>
void
Basis<T,Dual,R,CDF>::setLevel(int j) const
{
    _j = j;
}

template <typename T>
const BasisFunction<T,Dual,R,CDF> &
Basis<T,Dual,R,CDF>::generator(XType xtype) const
{
    if (xtype==XBSpline) {
        return mra.phi_;
    } else {
        return psi_;
    }
}

} // namespace lawa

#endif // LAWA_CONSTRUCTIONS_REALLINE_DUAL_BASIS_TCC
