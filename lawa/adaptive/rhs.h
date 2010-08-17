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


#ifndef ADAPTIVE_RHS_H
#define ADAPTIVE_RHS_H 1

#include <lawa/adaptive/index.h>
#include <lawa/adaptive/indexset.h>
#include <lawa/adaptive/coefficients.h>
#include <lawa/adaptive/referencesolutions.h>

namespace lawa {

template <typename T, typename Index, typename Basis, typename Preconditioner>
class RHS
{

};


template <typename T, typename Basis, typename Preconditioner>
class RHS<T,Index1D,Basis,Preconditioner>
{
	typedef typename Basis::BSplineType PrimalSpline;
	typedef typename Basis::WaveletType PrimalWavelet;

public:
	RHS(const Basis &basis, const Preconditioner &P, T (*f)(T), const DenseVector<Array<T> > &_singularPoints, const GeMatrix<flens::FullStorage<T, cxxblas::ColMajor> > &deltas);

	Coefficients<Lexicographical,T,Index1D>
	operator()(const IndexSet<Index1D> &Lambda);

	Coefficients<Lexicographical,T,Index1D>
	operator()(T tol);

private:
	const Basis &basis;
	const Preconditioner &P;
	PrimalSpline phi;
	PrimalWavelet psi;
	Function<T> f;
	const GeMatrix<flens::FullStorage<T, cxxblas::ColMajor> >	&deltas;
	Integral<T, Trapezoidal, PrimalSpline,  Function<T> > integral_sff;
	Integral<T, Trapezoidal, PrimalWavelet, Function<T> > integral_wf;
	Coefficients<Lexicographical,T,Index1D> rhs;
	Coefficients<AbsoluteValue,T,Index1D>   rhs_abs;


};


} // namespace lawa

#include <lawa/adaptive/rhs.tcc>

#endif // ADAPTIVE_RHS_H
