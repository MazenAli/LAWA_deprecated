/*
  LAWA - Library for Adaptive Wavelet Applications.
  Copyright (C) 2008,2009  Mario Rometsch, Alexander Stippler.

  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
*/

#ifndef APPLICATIONS_NEWEVALSCHEME_LOCSINGLESCALETRANSFORMS_H
#define APPLICATIONS_NEWEVALSCHEME_LOCSINGLESCALETRANSFORMS_H 1

#include <lawa/flensforlawa.h>
#include <lawa/constructions/basis.h>
#include <applications/new_eval_scheme/localrefinement.h>
#include <lawa/methods/adaptive/datastructures/datastructures.h>

namespace lawa {

template <typename PrimalBasis>
void
constructRandomGradedTree(const PrimalBasis &basis, int J, IndexSet<Index1D> &LambdaTree);

template <typename T, typename PrimalBasis>
void
computeLocalReconstruction(const Coefficients<Lexicographical,T,Index1D> &u_multi_j,
                           const PrimalBasis &basis, int j,
                           Coefficients<Lexicographical,T,Index1D> &u_multi_jP1);

template <typename T, typename PrimalBasis>
void
computeLocalReconstruction(const IndexSet<Index1D> &Lambda_multi_j,
                           const PrimalBasis &basis, int j,
                           IndexSet<Index1D> &Lambda_loc_single_jP1);

template <typename T, typename DualBasis>
void
computeLocalReconstruction_(const IndexSet<Index1D> &Lambda_multi_j,
                            const DualBasis &dual_basis, int j,
                            IndexSet<Index1D> &Lambda_loc_single_jP1);

template <typename T, typename PrimalBasis>
void
computeMultiToLocallySingleRepr(const PrimalBasis &basis,
                                const Coefficients<Lexicographical,T,Index1D> &u_multi,
                                Coefficients<Lexicographical,T,Index1D> &u_loc_single);

template <typename T, typename DualBasis>
void
computeLocalDecomposition(const Coefficients<Lexicographical,T,Index1D> &u_loc_single_j,
                          const DualBasis &dual_basis, int j,
                          const IndexSet<Index1D> &LambdaTree,
                          Coefficients<Lexicographical,T,Index1D> &u_loc_single_jM1,
                          Coefficients<Lexicographical,T,Index1D> &u_multi);

template <typename T, typename PrimalBasis>
void
computeLocalDecomposition_(const Coefficients<Lexicographical,T,Index1D> &u_loc_single_j,
                           const PrimalBasis &basis, int j,
                           const IndexSet<Index1D> &LambdaTree,
                           Coefficients<Lexicographical,T,Index1D> &u_loc_single_jM1,
                           Coefficients<Lexicographical,T,Index1D> &u_multi);

template <typename T, typename DualBasis>
void
computeLocallySingleToMultiRepr(const DualBasis &dual_basis,
                                const Coefficients<Lexicographical,T,Index1D> &u_loc_single,
                                const IndexSet<Index1D> &LambdaTree,
                                Coefficients<Lexicographical,T,Index1D> &u_multi);

template <typename T, typename PrimalBasis>
void
neighborhood(const PrimalBasis &basis, const Index1D &index, T c, IndexSet<Index1D> &ret,
             int shift=0);

template <typename T, typename PrimalBasis>
void
plot(const PrimalBasis &basis, const Coefficients<Lexicographical,T,Index1D> &u,
     std::stringstream &filename);


}   // namespace lawa

#include <applications/new_eval_scheme/loc_single_scale_transforms.tcc>

#endif // APPLICATIONS_NEWEVALSCHEME_LOCSINGLESCALETRANSFORMS_H
