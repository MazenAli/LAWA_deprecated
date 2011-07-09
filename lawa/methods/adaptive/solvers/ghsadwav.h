/*
  LAWA - Library for Adaptive Wavelet Applications.
  Copyright (C) 2008,2009  Sebastian Kestler, Mario Rometsch, Kristina Steih, Alexander Stippler.

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

#ifndef  LAWA_METHODS_ADAPTIVE_SOLVERS_GHSADWAV_H
#define  LAWA_METHODS_ADAPTIVE_SOLVERS_GHSADWAV_H 1

#include <lawa/methods/adaptive/datastructures/datastructures.h>
#include <lawa/methods/adaptive/algorithms/algorithms.h>

namespace lawa {

template <typename T, typename Index, typename AdaptiveOperator, typename RHS>
struct GHS_ADWAV {

        typedef typename IndexSet<Index>::const_iterator                       const_set_it;
        typedef typename Coefficients<Lexicographical,T,Index>::const_iterator const_coeff_it;
        typedef typename Coefficients<AbsoluteValue,T,Index>::const_iterator   const_coeff_abs_it;
        typedef typename Coefficients<Lexicographical,T,Index>::value_type     val_type;

        GHS_ADWAV(AdaptiveOperator &_A, RHS &_F);

        Coefficients<Lexicographical,T,Index>
        SOLVE(T nuM1, T _eps, int NumOfIterations=100, T H1norm=0.);

        IndexSet<Index>
        GROW(const Coefficients<Lexicographical,T,Index> &w, T nu_bar, T &nu);

        Coefficients<Lexicographical,T,Index>
        GALSOLVE(const IndexSet<Index> &Lambda, const Coefficients<Lexicographical,T,Index> &g,
                const Coefficients<Lexicographical,T,Index> &w, T delta, T tol);

        AdaptiveOperator &A;
        RHS              &F;
        T                cA, CA, kappa;
        T                alpha, omega, gamma, theta;
        T                eps;

        //        std::vector<Coefficients<Lexicographical,T,Index> > solutions;
        std::vector<T>                                      residuals;
        std::vector<T>                                      times;
        std::vector<int>                                    linsolve_iterations;

};


}    //namespace lawa

#include <lawa/methods/adaptive/solvers/ghsadwav.tcc>

#endif    // LAWA_METHODS_ADAPTIVE_SOLVERS_GHSADWAV_H

