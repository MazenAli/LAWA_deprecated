#ifndef LAWA_CONSTRUCTIONS_BASISFUNCTION_TCC
#define LAWA_CONSTRUCTIONS_BASISFUNCTION_TCC 1

#include <cassert>

namespace lawa {

template <typename T, FunctionSide Side, DomainType Domain, Construction Cons>
BasisFunction<T,Side,Domain,Cons>::~BasisFunction()
{
}

template <typename T, FunctionSide Side, DomainType Domain, Construction Cons>
T
BasisFunction<T,Side,Domain,Cons>::operator()(T /*x*/, int /*j*/, Integer /*k*/, 
                                              unsigned short /*deriv*/) const
{
    assert(0);
    return 0.;
}

template <typename T, FunctionSide Side, DomainType Domain, Construction Cons>
Support<T>
BasisFunction<T,Side,Domain,Cons>::support(int /*j*/, Integer /*k*/) const 
{
    assert(0);
    return Support<T>();
}

template <typename T, FunctionSide Side, DomainType Domain, Construction Cons>
DenseVector<Array<T> >
BasisFunction<T,Side,Domain,Cons>::singularSupport(int /*j*/, Integer /*k*/) const 
{
    assert(0);
    return DenseVector<Array<T> >(); 
}

template <typename T, FunctionSide Side, DomainType Domain, Construction Cons>
T
BasisFunction<T,Side,Domain,Cons>::tic(int /*j*/) const
{
    assert(0);
    return 0.;
}

} // namespace lawa

#endif // LAWA_CONSTRUCTIONS_BASISFUNCTION_H
