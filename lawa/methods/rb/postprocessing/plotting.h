#ifndef LAWA_METHODS_RB_POSTPROCESSING_PLOTTING_H
#define LAWA_METHODS_RB_POSTPROCESSING_PLOTTING_H

namespace lawa {

template <typename T, typename Basis>
void
saveCoeffVector2D(const Coefficients<Lexicographical,T,Index2D> &coeff, const Basis &basis2d, const char* filename);

template <typename T>
void
readCoeffVector2D(Coefficients<Lexicographical,T,Index2D>& coeff, const char* filename, bool append = false);


} // namespace lawa

#include <lawa/methods/rb/postprocessing/plotting.tcc>

#endif // LAWA_METHODS_RB_POSTPROCESSING_PLOTTING_H