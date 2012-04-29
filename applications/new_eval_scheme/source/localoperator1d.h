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

#ifndef APPLICATIONS_NEWEVALSCHEME_LOCALOPERATOR1D_H
#define APPLICATIONS_NEWEVALSCHEME_LOCALOPERATOR1D_H 1

#define DERIV 0

#include <lawa/flensforlawa.h>
#include <lawa/constructions/basis.h>
#include <lawa/methods/adaptive/algorithms/localrefinement2.h>

namespace lawa {

template <typename TestBasis, typename TrialBasis, typename BilinearForm>
class LocalOperator1D {

    public:
        typedef typename TrialBasis::T T;
        typedef TrialBasis                           TrialWaveletBasis;
        typedef TestBasis                            TestWaveletBasis;
        typedef typename TrialBasis::RefinementBasis TrialRefinementBasis;
        typedef typename TestBasis::RefinementBasis  TestRefinementBasis;

        typedef typename TreeCoefficients1D<T>::const_by_level_it                  const_by_level_it;
        typedef typename TreeCoefficients1D<T>::by_level_it                        by_level_it;

        LocalOperator1D(const TestBasis &_testBasis, const TrialBasis &_trialBasis,
                        BilinearForm &_Bil);

        const TestBasis                   &testBasis;
        const TrialBasis                  &trialBasis;
        const TestRefinementBasis         &testRefinementBasis;
        const TrialRefinementBasis        &trialRefinementBasis;
        BilinearForm                      &Bil;
        LocalRefinement2<TestBasis>       testLocalRefine;
        LocalRefinement2<TrialBasis>      trialLocalRefine;
        int                               testRefinementLevelOffset;
        int                               trialRefinementLevelOffset;

        void
        eval(const TreeCoefficients1D<T> &PsiLambdaHat, TreeCoefficients1D<T> &PsiLambdaCheck,
             const char* mode);

    private:
        void
        _evalA(int l, CoefficientsByLevel<T> &d, const TreeCoefficients1D<T> &c,
               CoefficientsByLevel<T> &PhiPiCheck, TreeCoefficients1D<T> &PsiLambdaCheck);

        void
        _evalA_nonRecursive(CoefficientsByLevel<T> &d, const TreeCoefficients1D<T> &c,
                            CoefficientsByLevel<T> &PhiPiCheck, TreeCoefficients1D<T> &PsiLambdaCheck);

        void
        _evalU(int l, CoefficientsByLevel<T> &d, const TreeCoefficients1D<T> &c,
               CoefficientsByLevel<T> &PhiPiCheck, TreeCoefficients1D<T> &PsiLambdaCheck);

        void
        _evalL(int l, CoefficientsByLevel<T> &d, const TreeCoefficients1D<T> &c,
               TreeCoefficients1D<T> &PsiLambdaCheck);

        void
        _splitPhiPi(int l, const CoefficientsByLevel<T> &c_l, CoefficientsByLevel<T> &PhiPiCheck1,
                    CoefficientsByLevel<T> &PhiPiCheck2) const;

        void
        _splitd(int l, const CoefficientsByLevel<T> &PsiLambdaCheck_l,
                CoefficientsByLevel<T> &d1, CoefficientsByLevel<T> &d2) const;

        void
        _applyBilinearForm(int l, const CoefficientsByLevel<T> &d,
                           CoefficientsByLevel<T> &PhiPiCheck);

};

}   // namespace lawa

#include <applications/new_eval_scheme/source/localoperator1d.tcc>

#endif // APPLICATIONS_NEWEVALSCHEME_LOCALOPERATOR1D_H
