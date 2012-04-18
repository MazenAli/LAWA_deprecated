#ifndef LAWA_CONSTRUCTIONS_INTERVAL_MULTI_REFINEMENTBASIS_TCC
#define LAWA_CONSTRUCTIONS_INTERVAL_MULTI_REFINEMENTBASIS_TCC 1

#include <cassert>
#include <lawa/math/linspace.h>
#include <lawa/math/pow2.h>
#include <lawa/constructions/interval/multi/_linear_evaluator.h>
#include <lawa/constructions/interval/multi/_quadratic_evaluator.h>
#include <lawa/constructions/interval/multi/_cubic_evaluator.h>

namespace lawa {

template <typename T>
Basis<T,Orthogonal,Interval,MultiRefinement>::Basis(int _d, int j)
    : mra(_d, j), d(_d), j0(mra.j0), _j(j0)
{
    assert(d>=2);
    setLevel(_j);

    this->enforceBoundaryCondition<DirichletBC>();
}

template <typename T>
Basis<T,Orthogonal,Interval,MultiRefinement>::~Basis()
{

}

template <typename T>
int
Basis<T,Orthogonal,Interval,MultiRefinement>::level() const
{
    return _j;
}

template <typename T>
void
Basis<T,Orthogonal,Interval,MultiRefinement>::setLevel(int j) const
{
    assert(j>=j0);
    _j = j;
}

template <typename T>
template <BoundaryCondition BC>
void
Basis<T,Orthogonal,Interval,MultiRefinement>::enforceBoundaryCondition()
{
    assert(BC==DirichletBC);
    mra.enforceBoundaryCondition<BC>();
}

template <typename T>
const BasisFunction<T,Orthogonal,Interval,MultiRefinement> &
Basis<T,Orthogonal,Interval,MultiRefinement>::generator(XType xtype) const
{
    if (xtype==XBSpline) {
        return mra.phi;
    } else {
        std::cerr << "BasisFunction<T,Orthogonal,Interval,MultiRefinement> has no wavelet member."
                  << std::endl;
        exit(1);
        return mra.phi;
    }
}

template <typename T>
template <typename SecondRefinementBasis>
void
Basis<T,Orthogonal,Interval,MultiRefinement>::
getBSplineNeighborsForBSpline(int j_bspline1, long k_bspline1,
                              const SecondRefinementBasis &secondrefinementbasis,
                              int &j_bspline2, long &k_bspline2_first, long &k_bspline2_last) const
{
    ct_assert(SecondRefinementBasis::Side==Orthogonal and SecondRefinementBasis::Domain==Interval
              and SecondRefinementBasis::Cons==MultiRefinement);
    //if (flens::IsSame<Basis<T,Orthogonal,Interval,MultiRefinement>, SecondRefinementBasis>::value)
    j_bspline2 = j_bspline1;
    Support<T> supp = mra.phi.support(j_bspline1,k_bspline1);
    T a = supp.l1, b = supp.l2;
    if (a==0.L) {
        k_bspline2_first = 0;       // In this case, range always starts with 0
        k_bspline2_last  = std::min(k_bspline1 + d, (long)mra.rangeI(j_bspline2).lastIndex());
        return;
    }
    if (b<1.L) {
        k_bspline2_first = std::max(k_bspline1 - d + 1, 0L);
        k_bspline2_last  = std::min(k_bspline1 + d - 1, (long)mra.rangeI(j_bspline2).lastIndex());
        return;
    }
    k_bspline2_first = std::max(0L,k_bspline1 - d);
    k_bspline2_last  = (long)mra.rangeI(j_bspline2).lastIndex();

    return;
}

template <typename T>
template <typename SecondBasis>
void
Basis<T,Orthogonal,Interval,MultiRefinement>
::getWaveletNeighborsForBSpline(int j_bspline, long k_bspline, const SecondBasis &secondbasis,
                                int &j_wavelet, long &k_wavelet_first, long &k_wavelet_last) const
{
    ct_assert(SecondBasis::Side==Orthogonal and SecondBasis::Domain==Interval
              and SecondBasis::Cons==Multi);

    j_wavelet = j_bspline - secondbasis._addRefinementLevel + 1;
    Support<T> supp = mra.phi.support(j_bspline,k_bspline);
    T a = supp.l1, b = supp.l2;

    if (a==0.L) {
        k_wavelet_first = 1;
        k_wavelet_last =  secondbasis._numLeftParts + secondbasis._numInnerParts;
        k_wavelet_last  = std::min(k_wavelet_last, (long)secondbasis.rangeJR(j_wavelet).lastIndex());
        return;
    }
    if (0<a && b<1.L) {
        long k_tilde = (long)std::floor(pow2i<T>(j_wavelet)*a);
        k_tilde += 1;
        T tmp = k_tilde*pow2i<T>(-j_wavelet);
        if (a<tmp && tmp<b) {
            k_wavelet_first  = (k_tilde-2)*secondbasis._numInnerParts;
            k_wavelet_last   = (k_tilde+2)*secondbasis._numInnerParts;
        }
        else {
            k_wavelet_first  = (k_tilde-2)*secondbasis._numInnerParts;
            k_wavelet_last   = (k_tilde+1)*secondbasis._numInnerParts;
        }
        k_wavelet_first = std::max((long)secondbasis.rangeJL(j_wavelet).firstIndex(),k_wavelet_first);
        k_wavelet_last  = std::min((long)secondbasis.rangeJR(j_wavelet).lastIndex(), k_wavelet_last);
        return;
    }
    k_wavelet_last   = secondbasis.rangeJ(j_wavelet).lastIndex();
    k_wavelet_first  = k_wavelet_last - (secondbasis._numRightParts + secondbasis._numInnerParts) + 1;
    k_wavelet_first  = std::max(1L, k_wavelet_first);
    return;
}

} // namespace lawa

#endif // LAWA_CONSTRUCTIONS_INTERVAL_MULTI_REFINEMENTBASIS_TCC
