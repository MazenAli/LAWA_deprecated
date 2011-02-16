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

#ifndef LAWA_ADAPTIVE_ALGORITHMS_GHSADWAV1D_H
#define LAWA_ADAPTIVE_ALGORITHMS_GHSADWAV1D_H 1

#include <lawa/adaptive/coefficients.h>
#include <lawa/adaptive/algorithms/apply1d.h>

namespace lawa {

template <typename T, typename Index, typename Basis1D, typename APPLY1D, typename RHS>
class GHS_ADWAV1D {

	typedef typename IndexSet<Index>::iterator 					    		set_it;
    typedef typename IndexSet<Index>::const_iterator 					    const_set_it;
    typedef typename Coefficients<Lexicographical,T,Index>::const_iterator  const_coeff_it;
    typedef typename Coefficients<AbsoluteValue,T,Index1D >::const_iterator const_coeff_abs_it;
    typedef typename Coefficients<Lexicographical,T,Index>::value_type      val_type;


public:

	GHS_ADWAV1D(const Basis1D &_basis, APPLY1D &_Apply, RHS &_F);

	Coefficients<Lexicographical,T,Index>
	SOLVE(T nuM1, T _eps, int NumOfIterations=100);

    std::vector<Coefficients<Lexicographical,T,Index> > solutions;
    std::vector<T>               residuals;
    std::vector<T>               times;


private:
	const Basis1D &basis;
	APPLY1D &Apply;
	RHS &F;
	T cA, CA, kappa;
	T alpha, omega, gamma, theta;
	T eps;

	IndexSet<Index>
	GROW(const Coefficients<Lexicographical,T,Index> &w, T nu_bar, T &nu);

	Coefficients<Lexicographical,T,Index>
	GALSOLVE(const IndexSet<Index> &Lambda, const Coefficients<Lexicographical,T,Index> &g,
			 const Coefficients<Lexicographical,T,Index> &w, T delta, T tol);
};


}	//namespace lawa

#include <lawa/adaptive/algorithms/ghsadwav1d.tcc>

#endif	//LAWA_ADAPTIVE_ALGORITHMS_GHSADWAV1D_H
