#ifndef LAWA_CONSTRUCTIONS_INTERVAL_DIJKEMA_DUAL_BSPLINE_TCC
#define LAWA_CONSTRUCTIONS_INTERVAL_DIJKEMA_DUAL_BSPLINE_TCC 1

#include <cassert>

namespace lawa {

template <typename T>
BSpline<T,Dual,Interval,Dijkema>::BSpline(const MRA<T,Dual,Interval,Dijkema> &_mra_)
    : mra_(_mra_)
{
}

template <typename T>
T
BSpline<T,Dual,Interval,Dijkema>::operator()(T x, int j, Integer k, unsigned short deriv) const
{
    assert(j>=mra_.j0);
    assert(k>=mra_.rangeI_(j).firstIndex());
    assert(k<=mra_.rangeI_(j).lastIndex());
    assert(deriv==0);

    if (k<=mra_.rangeI_L().lastIndex()) {
        T value = 0.0;
        Integer l = -mra_.phi_R.a_.lastIndex()+1;
        for ( Integer  r=mra_.R_Left.firstRow(); r<=mra_.R_Left.lastRow(); ++r, ++l) {
            value += mra_.R_Left(r,k) * mra_.phi_R(x,j,l);
        }
//        assert(l==-mra_.phi_R.a_.firstIndex());
        return value;
    }

    if (k>=mra_.rangeI_R(j).firstIndex()) {
        k -= mra_.rangeI_R(j).firstIndex()-1;
        T value = 0.0;
        Integer l = pow2i< Integer >(j)-mra_.phi_R.a_.lastIndex()+1;
        for ( Integer  r=mra_.R_Right.firstRow(); r<=mra_.R_Right.lastRow(); ++r, ++l) {
            value += mra_.R_Right(r,k) * mra_.phi_R(x,j,l);
        }
        return value;
    }

    return mra_.phi_R(x,j,k-(mra_.d + mra_.d_ - 1)-mra_.phi_R.l1_);
}

template <typename T>
Support<T>
BSpline<T,Dual,Interval,Dijkema>::support(int j, Integer k) const
{
    if (k<=mra_.rangeI_L().lastIndex()) {
        return Support<T>(0.,pow2i<T>(-j)*(mra_.phi_R.a_.length()-2));
    }
    if (k>=mra_.rangeI_R(j).firstIndex()) {
        return Support<T>(1.-pow2i<T>(-j)*(mra_.phi_R.a_.length()-2), 1.);
    }
    return pow2i<T>(-j)*Support<T>(k-(mra_.d+mra_.d_-1),
                                   k-(mra_.d+mra_.d_-1)+mra_.phi_R.l2_-mra_.phi_R.l1_);
}

} // namespace lawa

#endif // LAWA_CONSTRUCTIONS_INTERVAL_DIJKEMA_DUAL_BSPLINE_TCC
