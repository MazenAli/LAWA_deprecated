#ifndef LAWA_METHODS_RB_POSTPROCESSING_PLOTTING_H
#define LAWA_METHODS_RB_POSTPROCESSING_PLOTTING_H

#include <lawa/methods/adaptive/datastructures/datastructures.h>

namespace lawa {

// Version NonPeriodic x NonPeriodic
template <typename T, typename Basis>
typename RestrictTo<SFINAE_Wrapper<!IsPeriodic<typename Basis::FirstBasisType>::value
					and !IsPeriodic<typename Basis::SecondBasisType>::value, T>::value, void>::Type
saveCoeffVector2D(const Coefficients<Lexicographical,T,Index2D> &coeff, const Basis &basis2d, const char* filename);

// Version Periodic x NonPeriodic
template <typename T, typename Basis>
typename RestrictTo<SFINAE_Wrapper<IsPeriodic<typename Basis::FirstBasisType>::value
					and !IsPeriodic<typename Basis::SecondBasisType>::value, T>::value, void>::Type
saveCoeffVector2D(const Coefficients<Lexicographical,T,Index2D> &coeff, const Basis &basis2d, const char* filename);

template <typename T>
void
readCoeffVector2D(Coefficients<Lexicographical,T,Index2D>& coeff, const char* filename, bool append = false);

template <typename T, typename Basis>
void
saveIndexSet2D(const IndexSet<Index2D> &indexset, const Basis &basis2d, const char* filename);
  
template <typename T>
void
readIndexSet2D(IndexSet<Index2D>& indexset, const char* filename, bool append = false);

} // namespace lawa

#include <lawa/methods/rb/postprocessing/plotting.tcc>

#endif // LAWA_METHODS_RB_POSTPROCESSING_PLOTTING_H
