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


#ifndef ADAPTIVE_BILINEARFORM_H
#define ADAPTIVE_BILINEARFORM_H 1

#include <adaptive/index.h>
#include <adaptive/bilinearform.h>
#include <lawa/lawa.h>

namespace lawa {

/* Minimal structure of a BilinearForm

template <typename T, typename Index, typename Basis>
class BilinearForm
{
public:
	BilinearForm(const BilinearForm<T,Index,Basis> &_a);

	T
	operator()(const Index &row_index, const Index &col_index);
};

*/

template <typename T, typename Basis>
class HelmholtzOperator1d
{
private:
	const Basis& basis;
	const T c;	//reaction term

	typedef typename Basis::BSplineType PrimalSpline;
	typedef typename Basis::WaveletType PrimalWavelet;

	PrimalSpline  phi, d_phi;
	PrimalWavelet psi, d_psi;

    Integral<T, Gauss, PrimalSpline, PrimalSpline>   integral_sfsf, dd_integral_sfsf;
    Integral<T, Gauss, PrimalSpline, PrimalWavelet>  integral_sfw,  dd_integral_sfw;
    Integral<T, Gauss, PrimalWavelet, PrimalSpline>  integral_wsf,  dd_integral_wsf;
    Integral<T, Gauss, PrimalWavelet, PrimalWavelet> integral_ww,   dd_integral_ww;

public:
	HelmholtzOperator1d(const Basis &basis, const T &c);
	HelmholtzOperator1d(const HelmholtzOperator1d<T,Basis> &a);

	T getc() const;
	const Basis& getBasis() const;

	T
	operator()(const Index1d &row_index, const Index1d &col_index) const;
};


}	// namespace lawa

#include <adaptive/bilinearform.tcc>

#endif // ADAPTIVE_BILINEARFORM_H
