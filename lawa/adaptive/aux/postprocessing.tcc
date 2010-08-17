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

namespace lawa {

template <typename T, typename Index, typename MA, typename RHS>
T
estimateError_Au_M_f(MA &A, RHS &F, const Coefficients<Lexicographical,T,Index> & u,
					 const IndexSet<Index> &LambdaCol)
{
	Coefficients<Lexicographical,T,Index> Au(u.d,u.d_), f(u.d,u.d_), res(u.d,u.d_);
	Au = mv(LambdaCol,A,u);
	f  = F(LambdaCol);
	res = Au-f;
	return res.norm(2.);
}

template <typename T, typename Index, typename MA, typename RHS>
T
estimateError_H_energy(MA &A_H, RHS &F_H, const Coefficients<Lexicographical,T,Index> & u,
					   T HNormOfExactSolution)
{
	Coefficients<Lexicographical,T,Index> Au(u.d,u.d_), f_H(u.d,u.d_);
	Au = mv(supp(u),A_H,u);
	T uAu = u*Au;
	f_H   = F_H(supp(u));
	T fu  = f_H*u;

	return std::sqrt(fabs(std::pow(HNormOfExactSolution,2)- 2*fu + uAu));

}

} // namespace lawa
