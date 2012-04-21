#ifndef LAWA_CONSTRUCTIONS_SPARSEMULTI_INTERVAL__SPARSEMULTI_SCALING_EVALUATOR_H
#define LAWA_CONSTRUCTIONS_SPARSEMULTI_INTERVAL__SPARSEMULTI_SCALING_EVALUATOR_H 1

namespace lawa {

//--- cubic evaluators -----------------------------------------------------

template <typename T>
    T
    _cubic_sparsemulti_scaling_left_evaluator0(T x, unsigned short deriv);

template <typename T>
    T
    _cubic_sparsemulti_scaling_inner_evaluator0(T x, unsigned short deriv);

template <typename T>
    T
    _cubic_sparsemulti_scaling_inner_evaluator1(T x, unsigned short deriv);

template <typename T>
    T
    _cubic_sparsemulti_scaling_right_evaluator0(T x, unsigned short deriv);

} // namespace lawa

#include <lawa/constructions/sparsemulti/interval/_sparsemulti_scaling_evaluator.tcc>

#endif // LAWA_CONSTRUCTIONS_SPARSEMULTI_INTERVAL__SPARSEMULTI_SCALING_EVALUATOR_H
