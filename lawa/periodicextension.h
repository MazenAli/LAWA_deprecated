#ifndef LAWA_PERIODICEXTENSION_H
#define LAWA_PERIODICEXTENSION_H 1

#include <lawa/flensforlawa.h>

namespace lawa {

using namespace flens;

template <typename T>
struct PeriodicExtension
{
    PeriodicExtension(Function<T> _F);

    T
    operator()(T x) const;

    DenseVector<Array<T> > singularPoints;
    T (*f)(T);
};

} // namespace lawa

#include <lawa/periodicextension.tcc>

#endif // LAWA_PERIODICEXTENSION_H