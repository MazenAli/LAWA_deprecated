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

} // namespace lawa

#endif // LAWA_CONSTRUCTIONS_INTERVAL_MULTI_REFINEMENTBASIS_TCC
