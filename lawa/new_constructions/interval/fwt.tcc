#ifndef LAWA_CONSTRUCTIONS_INTERVAL_FWT_TCC
#define LAWA_CONSTRUCTIONS_INTERVAL_FWT_TCC 1

#include <cassert>

namespace lawa {

template <typename X, typename Y, Construction Cons>
void
decompose(const DenseVector<X> &x, 
          const Basis<typename X::ElementType,Dual,Interval,Cons> &basis_, int j,
          DenseVector<Y> &y)
{
    assert(j>=basis_.j0);
    assert(x.range()==basis_.mra_.rangeI_(j+1));
    typedef typename X::ElementType T;

    basis_.setLevel(j);
    DenseVector<Array<T> > sf, w;
    sf = transpose(basis_.mra_.M0_) * x;
    w  = transpose(basis_.M1_) * x;
    concatenate(sf,w,y,y.firstIndex());
}

template <typename X, typename Y, Construction Cons>
void
reconstruct(const DenseVector<X> &x, 
            const Basis<typename X::ElementType,Primal,Interval,Cons> &basis, int j,
            DenseVector<Y> &y)
{
    typedef typename X::ElementType T;
    assert(j>=basis.j0);
    assert(x.range()==basis.mra.rangeI(j+1));
    basis.setLevel(j);
    y  = basis.mra.M0 * x(basis.mra.rangeI(j));
    y += basis.M1*x(_(basis.mra.rangeI(j).lastIndex()+1,basis.mra.rangeI(j+1).lastIndex()));
}

template <typename X, typename Y, Construction Cons>
void
fwt(const DenseVector<X> &x, 
    const Basis<typename X::ElementType,Dual,Interval,Cons> &basis_, int j,
    DenseVector<Y> &y)
{
    assert(j>=basis_.j0);
    assert(x.range()==basis_.mra_.rangeI_(j+1));
    
    y = x;
    for (int l=j; l>=basis_.j0; --l) {
        typename DenseVector<Y>::View yview = y(basis_.mra_.rangeI_(l+1));
        decompose(yview, basis_, l, yview);
    }
}

template <typename X, typename Y, Construction Cons>
void
ifwt(const DenseVector<X> &x, 
     const Basis<typename X::ElementType,Primal,Interval,Cons> &basis, int j,
     DenseVector<Y> &y)
{
    assert(j>=basis.j0);
    assert(x.range()==basis.mra.rangeI(j+1));
    
    y = x;
    for (int l=basis.j0; l<=j; ++l) {
        typename DenseVector<Y>::View yview = y(basis.mra.rangeI(l+1));
        DenseVector<Array<typename X::ElementType> > z = y(basis.mra.rangeI(l+1));
        reconstruct(z, basis, l, yview);
    }
}

} // namespace lawa

#endif // LAWA_CONSTRUCTIONS_INTERVAL_FWT_H
