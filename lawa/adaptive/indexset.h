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

#ifndef LAWA_ADAPTIVE_INDEXSET_H_
#define LAWA_ADAPTIVE_INDEXSET_H_ 1

#include <set>
#include <lawa/adaptive/index.h>
#include <lawa/basis.h>
#include <lawa/bspline.h>
#include <lawa/enum.h>
#include <lawa/mra.h>
#include <lawa/wavelet.h>

namespace lawa {

template <typename Index>
struct IndexSet : std::set<Index, lt<Lexicographical, Index > >
{
    IndexSet(int d, int d_);

    IndexSet<Index>&
    operator= (const IndexSet<Index> &_set);

    IndexSet<Index>
    operator+ (const IndexSet<Index> &_set) const;

    const int d,d_;
};

template <typename Index>
std::ostream& operator<< (std::ostream &s, const IndexSet<Index> &i);

void
getMinAndMaxLevel(const IndexSet<Index1D> &Lambda, int &jmin, int &jmax);

IndexSet<Index1D>
extractSpaceIndices(const IndexSet<Index2D> &Lambda);

//Security zone for an 1d-index following Urban:2009, p.235 and KU:2010.
template <typename T, DomainType Domain, Construction Cons>
IndexSet<Index1D>
C(const IndexSet<Index1D> &Lambda, T c, const Basis<T,Primal,Domain,Cons> &basis);

template <typename T, DomainType Domain, Construction Cons>
IndexSet<Index1D>
C(const Index1D &lambda, T c, const Basis<T,Primal,Domain,Cons> &basis);

template <typename T>
IndexSet<Index1D>
C_WO_XBSpline(const IndexSet<Index1D> &Lambda, T c, const Basis<T,Primal,R,CDF> &basis);

template <typename T>
IndexSet<Index1D>
C_WO_XBSpline(const Index1D &lambda, T c, const Basis<T,Primal,R,CDF> &basis);


template <typename T, DomainType Domain, Construction Cons>
void
C(const Index1D &lambda, T c, const MRA<T,Primal,Domain,Cons> &mra,
  const Basis<T,Primal,Domain,Cons> &basis, IndexSet<Index1D> &ret);

// Computation of a security zone for 2d-tensor basis
template <typename T, typename Basis2D>
IndexSet<Index2D>
C(const IndexSet<Index2D> &Lambda, T c, const Basis2D &basis);

template <typename T, typename Basis2D>
IndexSet<Index2D>
C_t(const IndexSet<Index2D> &Lambda, T c, const Basis2D &basis);

template <typename T, typename Basis3D>
IndexSet<Index3D>
C(const IndexSet<Index3D> &Lambda, T c, const Basis3D &basis);


// Realizations of lambdaTilde for different Basis
template <typename T, typename Basis>
IndexSet<Index1D>
lambdaTilde1d_PDE(const Index1D &lambda, const Basis &basis,
                  int s_tilde, int jmin, int jmax, bool update);


}   // namespace lawa


#include <lawa/adaptive/indexset.tcc>


#endif // LAWA_ADAPTIVE_INDEXSET_H_ */
