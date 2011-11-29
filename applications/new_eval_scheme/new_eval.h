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

#ifndef APPLICATIONS_NEWEVALSCHEME_NEWEVAL_H
#define APPLICATIONS_NEWEVALSCHEME_NEWEVAL_H 1

#include <lawa/flensforlawa.h>
#include <lawa/constructions/basis.h>
#include <lawa/methods/adaptive/datastructures/datastructures.h>
#include <applications/new_eval_scheme/loc_single_scale_transforms.h>

namespace lawa {

template <typename T, typename PrimalTestBasis, typename PrimalTrialBasis>
void
eval(int l, const PrimalTestBasis &test_basis, const PrimalTrialBasis &trial_basis,
            std::map<int, Coefficients<Lexicographical,T,Index1D> > &c,
            const Coefficients<Lexicographical,T,Index1D> &d_lM1,
            std::map<int, IndexSet<Index1D> > &Lhd,
            const IndexSet<Index1D> &SquareCap_lM1,
            Coefficients<Lexicographical,T,Index1D> &UpsilonVsV_lM1,
            Coefficients<Lexicographical,T,Index1D> &ThetaVsV_lM1);

}   // namespace lawa

#include <applications/new_eval_scheme/new_eval.tcc>

#endif  // APPLICATIONS_NEWEVALSCHEME_NEWEVAL_H
