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

#ifndef APPLICATIONS_NEWEVALSCHEME_LOCALOPERATOR_H
#define APPLICATIONS_NEWEVALSCHEME_LOCALOPERATOR_H 1

#define DERIV 0

#include <lawa/flensforlawa.h>
#include <lawa/constructions/basis.h>
#include <lawa/methods/adaptive/algorithms/localrefinement.h>

namespace lawa {

template <typename TestBasis, typename TrialBasis, typename BilinearForm, typename Preconditioner>
struct LocalOperator {

    typedef typename TrialBasis::T T;

    typedef IndexSet<Index1D>::const_iterator                                  const_set1d_it;
    typedef typename Coefficients<Lexicographical,T,Index1D>::const_iterator   const_coeff1d_it;
    typedef typename TreeCoefficients1D<T>::const_by_level_it                  const_by_level_it;
    typedef typename TreeCoefficients1D<T>::by_level_it                        by_level_it;

    LocalOperator(const TestBasis &_test_basis,   bool test_withDirichletBC,
                  const TrialBasis &_trial_basis, bool trial_withDirichletBC,
                  const int offset,
                  const BilinearForm &_Bil, const Preconditioner &_Prec);

    void
    scale_wrt_trialbasis(const TreeCoefficients1D<T> &x, TreeCoefficients1D<T> &y);

    void
    evalA(int l, const CoefficientsByLevel<T> &d, const TreeCoefficients1D<T> &c,
          CoefficientsByLevel<T> &PhiPiCheck_vs_v, TreeCoefficients1D<T> &PsiLambdaCheck_vs_v,
          bool pre_apply_prec=true);

    void
    evalU(int l, const CoefficientsByLevel<T> &d, const TreeCoefficients1D<T> &c,
          CoefficientsByLevel<T> &PhiPiCheck_vs_v, TreeCoefficients1D<T> &PsiLambdaCheck_vs_v,
          bool pre_apply_prec=true);

    void
    evalL(int l, const CoefficientsByLevel<T> &d, const TreeCoefficients1D<T> &c,
          TreeCoefficients1D<T> &PsiLambdaCheck_vs_v, bool pre_apply_prec=true);


    const TestBasis                   &test_basis;
    const TrialBasis                  &trial_basis;
    const BilinearForm                &Bil;
    const Preconditioner              &Prec;
    LocalRefinement<TestBasis>        test_localtransform;
    LocalRefinement<TrialBasis>       trial_localtransform;

    //Important for integration: we only into account indices within the range {k-offset,...,k+offset}.
    //This approach might fail if offset is not chosen w.r.t. the underlying bases!!
    const int                         offset;
};

}   // namespace lawa

#include <applications/new_eval_scheme/source/localoperator.tcc>

#endif // APPLICATIONS_NEWEVALSCHEME_LOCALOPERATOR_H
