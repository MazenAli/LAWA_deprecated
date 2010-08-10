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

#ifndef INDEXSET_H_
#define INDEXSET_H_ 1

#include <set>
#include <lawa/enum.h>
#include <adaptive/index.h>
#include <lawa/lawa.h>

namespace lawa {

template <typename Index>
struct IndexSet : std::set<Index, lt<Lexicographical,Index > >
{
    IndexSet<Index>&
    operator= (const IndexSet<Index> &_set);

    IndexSet<Index>
    operator+ (const IndexSet<Index> &_set) const;
};

template <typename Index>
std::ostream& operator<< (std::ostream &s, const IndexSet<Index> &i);

template <typename T, DomainType Domain, Construction Cons>
IndexSet<Index1d>
C(const Index1d &lambda, T c, const BSpline<T,Primal,Domain,Cons> &phi);

template <typename T, DomainType Domain, Construction Cons>
IndexSet<Index1d>
C(const Index1d &lambda, T c, const Wavelet<T,Primal,Domain,Cons> &psi);

template <typename T>
IndexSet<Index1d>
C_interval(const IndexSet<Index1d> &Lambda, T c);

template <typename T>
IndexSet<Index1d>
C_realline(const IndexSet<Index1d> &Lambda, T c);

template <typename T>
IndexSet<Index1d>
C_periodic(const IndexSet<Index1d> &Lambda, T c);

}   // namespace lawa


#include "indexset.tcc"


#endif /* INDEXSET_H_ */