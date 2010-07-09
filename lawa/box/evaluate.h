#ifndef LAWA_BOX_EVALUATE_H
#define LAWA_BOX_EVALUATE_H 1

namespace lawa {

template <typename X, typename Basis>
typename X::ElementType
evaluate(const Basis &basis, const int J_x, const int J_y, const flens::DenseVector<X> &coeffs,
         const typename X::ElementType x, const typename X::ElementType y, const int deriv_x,
         const int deriv_y);

} // namespace lawa

#include <lawa/box/evaluate.tcc>


#endif // LAWA_BOX_EVALUATE_H