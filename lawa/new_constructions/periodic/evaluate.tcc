#ifndef LAWA_CONSTRUCTIONS_PERIODIC_EVALUATE_TCC
#define LAWA_CONSTRUCTIONS_PERIODIC_EVALUATE_TCC 1

#include <lawa/aux/integer.h>

namespace lawa {
    

template <typename X>
typename X::ElementType
evaluate(const MRA<typename X::ElementType,Primal,Periodic,CDF> &mra, int j,
         const DenseVector<X> &coeffs, typename X::ElementType x, int deriv)
{
    typedef typename X::ElementType T;
    assert(j>=mra.j0);
    assert(coeffs.length()==mra.cardI(j));
    assert(x>=0.);
    assert(x<=1.);
    
    BSpline<T,Primal,Periodic,CDF> phi(mra.d);
    int offsetI = mra.rangeI(mra.j0).firstIndex()-coeffs.firstIndex();
    T ret = 0.0;
    for (int k=mra.rangeI(j).firstIndex(); k<=mra.rangeI(j).lastIndex(); ++k) {
        ret += coeffs(k-offsetI) * phi(x,j,k,deriv);
    }
    return ret;
}

template <typename X>
typename X::ElementType
evaluate(const Basis<typename X::ElementType,Primal,Periodic,CDF> &basis,
         int J, const DenseVector<X> &coeffs, typename X::ElementType x, 
         int deriv)
{
    typedef typename X::ElementType T;
    assert(J>=basis.j0);
    assert(coeffs.length()==basis.mra.cardI(J));
    assert(x>=0.);
    assert(x<=1.);

    const int j0 = basis.j0;
    basis.setLevel(j0);
    T ret = 0;
    int offsetJ = basis.rangeJ(j0).firstIndex()-coeffs.firstIndex();
    Range<Integer> range(coeffs.firstIndex(), coeffs.firstIndex() + basis.cardJ(j0) - 1);
    
    ret += evaluate(basis.mra,j0,coeffs(range),x,deriv);
    Wavelet<T,Primal,Periodic,CDF> psi(basis.d, basis.d_);
    for (int j=j0; j<=J-1; ++j) {
        for (int k=basis.rangeJ(j).firstIndex(); k<=basis.rangeJ(j).lastIndex(); ++k) {
            ret += coeffs(basis.mra.cardI(j) + k - offsetJ) * basis.psi(x,j,k,deriv);
        }
    }
    return ret;
}

} // namespace lawa

#endif // LAWA_CONSTRUCTIONS_PERIODIC_EVALUATE_TCC
