#ifndef LAWA_CONSTRUCTIONS_SPARSEMULTI_INTERVAL__SPARSEMULTI_REALLINE_SCALING_EVALUATOR_H
#define LAWA_CONSTRUCTIONS_SPARSEMULTI_INTERVAL__SPARSEMULTI_REALLINE_SCALING_EVALUATOR_H 1

#include <lawa/constructions/sparsemulti/interval/_sparsemulti_scaling_evaluator.h>

namespace lawa {

//--- cubic evaluators -----------------------------------------------------

template <typename T>
    T
    _cubic_realline_sparsemulti_scaling_inner_evaluator1(T x, unsigned short deriv);

template <typename T>
    T
    _cubic_realline_sparsemulti_scaling_right_evaluator0(T x, unsigned short deriv);

} // namespace lawa

#include <lawa/constructions/sparsemulti/realline/_sparsemulti_realline_scaling_evaluator.tcc>

#endif // LAWA_CONSTRUCTIONS_SPARSEMULTI_INTERVAL__REALLINE_SPARSEMULTI_SCALING_EVALUATOR_H
