/*
 * indexset_generation.h
 *
 *  Created on: 26.02.2013
 *      Author: ksteih
 */

#ifndef LAWA_METHODS_ADAPTIVE_ALGORITHMS_INDEXSET_GENERATION_H_
#define LAWA_METHODS_ADAPTIVE_ALGORITHMS_INDEXSET_GENERATION_H_

#include <lawa/methods/adaptive/datastructures/datastructures.h>

namespace lawa {

template <typename T, typename Basis2D>
void
getSparseGridIndexSet(const Basis2D &basis, IndexSet<Index2D> &Lambda, int j, int deltaL, T gamma = 0.);

} // namespace lawa

#include <lawa/methods/adaptive/algorithms/indexset_generation.tcc>

#endif /* LAWA_METHODS_ADAPTIVE_ALGORITHMS_INDEXSET_GENERATION_H_ */
