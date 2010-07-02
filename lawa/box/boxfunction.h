#ifndef LAWA_BOX_BOXFUNCTION_H
#define LAWA_BOX_BOXFUNCTION_H 1

#include <lawa/flensforlawa.h>
#include <vector>

namespace lawa{
    
using namespace flens;

template <typename T>
struct BoxFunction
{
    BoxFunction(T (*_f)(T,T), std::vector<DenseVector<Array<T> > > _singularPoints );
    
    T
    operator()(T x, T y) const;
    
    T (*f)(T,T);
    std::vector<DenseVector<Array<T> > > singularPoints;
};    
    
} // namespace lawa

#include <lawa/box/boxfunction.tcc>

#endif // LAWA_BOX_BOXFUNCTION_H