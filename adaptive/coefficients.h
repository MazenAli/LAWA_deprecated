/*
  LAWA - Library for Adaptive Wavelet Applications.
  Copyright (C) 2008,2009  Mario Rometsch, Kristina Steih, Alexander Stippler.

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

#ifndef COEFFICIENTS_H
#define COEFFICIENTS_H 1

#include <set>
#include <adaptive/index.h>
#include <adaptive/indexset.h>
#include <lawa/lawa.h>

namespace lawa {

template <SortingCriterion S, typename T, typename Index>
struct Coefficients
{
};

template <typename T, typename Index>
struct Coefficients<Lexicographical,T,Index> : std::map<Index,T,lt<Lexicographical,Index> >
{
	Coefficients();		//todo: where is this required?

	Coefficients(const int d, const int d_);

	Coefficients<Lexicographical,T,Index>&
	operator=(const Coefficients<Lexicographical,T,Index> &_coeff);

	Coefficients<Lexicographical,T,Index>&
	operator=(const Coefficients<AbsoluteValue,T,Index> &_coeff);

	Coefficients<Lexicographical,T,Index>
	operator-(const Coefficients<Lexicographical,T,Index> &_coeff) const;

	Coefficients<Lexicographical,T,Index>
	operator+(const Coefficients<Lexicographical,T,Index> &_coeff) const;

	T
	operator*(const Coefficients<Lexicographical,T,Index> &_coeff) const;

	T
	norm(T tau=2.0) const;

	const int d, d_;
};

template <typename T, typename Index>
std::ostream& operator<< (std::ostream &s, const Coefficients<Lexicographical,T,Index> &c);

template <typename T, typename Index>
Coefficients<Lexicographical,T,Index>
operator*(T alpha, const Coefficients<Lexicographical,T,Index> &_coeff);

template <typename T, typename Index>
IndexSet<Index>
supp(const Coefficients<Lexicographical,T,Index> &v);

template <typename T, typename Index>
void
FillWithZeros(const IndexSet<Index> &Lambda, Coefficients<Lexicographical,T,Index> &_coeff);


template <typename T, typename Index>
struct Coefficients<AbsoluteValue,T,Index> : std::multimap<T,Index,lt<AbsoluteValue,T> >
{
	Coefficients(const int d, const int d_);

	Coefficients<AbsoluteValue,T,Index>&
	operator=(const Coefficients<Lexicographical,T,Index> &_coeff);

	Coefficients<AbsoluteValue,T,Index>&
	operator=(const Coefficients<AbsoluteValue,T,Index> &_coeff);

	T
	norm(T tau=2.0) const;

	T
	wtauNorm(T tau) const;

	DenseVector<Array<T> >
	norm_sections() const;

	const int d, d_;
};

template <typename T, typename Index>
std::ostream& operator<< (std::ostream &s, const Coefficients<AbsoluteValue,T,Index> &c);


} // namespace lawa

#include <adaptive/coefficients.tcc>

#endif // COEFFICIENTS_H
