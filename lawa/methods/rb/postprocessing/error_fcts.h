#ifndef LAWA_METHODS_RB_POSTPROCESSING_ERROR_FCTS_H
#define LAWA_METHODS_RB_POSTPROCESSING_ERROR_FCTS_H

namespace lawa {

template<typename T, typename Basis2D, typename Prec>
T
L2_H1_error(Basis2D& basis2d, const Coefficients<Lexicographical,T,Index2D> &u, Prec& P,
            T (*u_ref)(T,T), T (*dx_u_ref)(T,T), T a_t, T b_t, int n_t, T a_x, T b_x, int n_x);

} // namespace lawa

#include <lawa/methods/rb/postprocessing/error_fcts.tcc>

#endif // LAWA_METHODS_RB_POSTPROCESSING_ERROR_FCTS_H
