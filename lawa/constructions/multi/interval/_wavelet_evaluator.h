#ifndef LAWA_CONSTRUCTIONS_MULTI_INTERVAL__WAVELET_EVALUATOR_H
#define LAWA_CONSTRUCTIONS_MULTI_INTERVAL__WAVELET_EVALUATOR_H 1

namespace lawa {

//--- linear evaluators --------------------------------------------------------

template <typename T>
    T
    _linear_wavelet_left_evaluator0(T x, unsigned short deriv);

template <typename T>
    T
    _linear_wavelet_left_evaluator1(T x, unsigned short deriv);

template <typename T>
    T
    _linear_wavelet_inner_evaluator0(T x, unsigned short deriv);

template <typename T>
    T
    _linear_wavelet_inner_evaluator1(T x, unsigned short deriv);

template <typename T>
    T
    _linear_wavelet_inner_evaluator2(T x, unsigned short deriv);

template <typename T>
    T
    _linear_wavelet_right_evaluator0(T x, unsigned short deriv);

template <typename T>
    T
    _linear_wavelet_right_evaluator1(T x, unsigned short deriv);

} // namespace lawa

#include <lawa/constructions/multi/interval/_wavelet_evaluator.tcc>

#endif // LAWA_CONSTRUCTIONS_MULTI_INTERVAL__WAVELET_EVALUATOR_H
