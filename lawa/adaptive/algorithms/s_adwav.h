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
#ifndef ADAPTIVE_ALGORITHMS_S_ADWAV_H
#define ADAPTIVE_ALGORITHMS_S_ADWAV_H 1

#include <iostream>
#include <vector>
#include <lawa/adaptive/indexset.h>
#include <lawa/adaptive/coefficients.h>
#include <lawa/adaptive/mapmatrix.h>
#include <lawa/adaptive/aux/postprocessing.h>

namespace lawa {

template <typename T, typename Index, typename Basis, typename MA, typename RHS>
class S_ADWAV {
public:
    S_ADWAV(const Basis &basis, MA &A, RHS &F, T contraction, T start_threshTol,
            T start_linTol, T start_resTol, int max_iters, T eps);

	void solve_cg(const IndexSet<Index> &Initial_Lambda);
	void solve_cg_with_error_on_the_fly(const IndexSet<Index> &Initial_Lambda, T H1norm);
	void solve_gmres(const IndexSet<Index> &Initial_Lambda);

    std::vector<Coefficients<Lexicographical,T,Index> > solutions;
    std::vector<T>               residuals;
    std::vector<T>               times;
    
private:
    const Basis &basis;
    MA &A;
    RHS &F;
    T contraction, threshTol, linTol, resTol;
    int NumOfIterations;
    T eps;

};

} //namespace lawa

#include <lawa/adaptive/algorithms/s_adwav.tcc>

#endif //ADAPTIVE_ALGORITHMS_S_ADWAV_H
