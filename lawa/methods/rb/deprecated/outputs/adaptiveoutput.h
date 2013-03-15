#ifndef LAWA_METHODS_RB_OUTPUTS_ADAPTIVEOUTPUT_H
#define LAWA_METHODS_RB_OUTPUTS_ADAPTIVEOUTPUT_H 1

#include <lawa/methods/adaptive/datastructures/adaptiverhs.h>

namespace lawa {

template<typename T, typename Index>
struct AdaptiveOutput : public AdaptiveRhs<T, Index> {

    virtual T
    operator()(const Coefficients<Lexicographical, T, Index>& coeffs_u) = 0;
    
};
    
} // namespace lawa

#endif // LAWA_METHODS_RB_OUTPUTS_ADAPTIVEOUTPUT_H
