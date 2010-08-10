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

#ifndef INDEXSET_TCC_
#define INDEXSET_TCC_

namespace lawa {


template <typename Index>
IndexSet<Index>&
IndexSet<Index>::operator=(const IndexSet<Index> &_set)
{
    typedef typename IndexSet<Index>::const_iterator const_it;
    erase(IndexSet<Index>::begin(), IndexSet<Index>::end());
    if (_set.size() > 0) {
        for (const_it lambda = _set.begin(); lambda != _set.end(); ++lambda) {
            insert(*lambda);
        }
    }
    return *this;
}

template <typename Index>
IndexSet<Index>
IndexSet<Index>::operator+(const IndexSet<Index> &_set) const
{
    typedef typename IndexSet<Index>::const_iterator const_it;
    IndexSet<Index> ret = *this;
    for (const_it lambda = _set.begin(); lambda != _set.end(); ++lambda) {
        ret.insert(*lambda);
    }
    return ret;
}

template <typename Index>
std::ostream& operator<< (std::ostream &s, const IndexSet<Index> &i)
{
    typedef typename IndexSet<Index>::const_iterator const_it;
    s << std::endl << "IndexSet<Index>:" << std::endl;
    if (i.size() > 0) {
        for (const_it lambda = i.begin(); lambda != i.end(); ++lambda) {
            s << "  [" << (*lambda) << "]" << std::endl;
        }
    }
    return s << std::endl;
}

template <typename T>
IndexSet<Index1d>
C(const Index1d &lambda, T c, const BSpline<T,Primal,R,CDF> &basis) {
	std::cout << "Security zone for realline BSpline!" << std::endl;
	IndexSet<Index1d> indexset;
	return indexset;
}

template <typename T>
IndexSet<Index1d>
C(const Index1d &lambda, T c, const Wavelet<T,Primal,R,CDF> &basis) {
	std::cout << "Security zone for realline wavelet!" << std::endl;
	IndexSet<Index1d> indexset;
	return indexset;
}

template <typename T>
IndexSet<Index1d>
C_interval(const IndexSet<Index1d> &Lambda, T c) {
	IndexSet<Index1d> indexset;
	return indexset;
}

template <typename T>
IndexSet<Index1d>
C_realline(const IndexSet<Index1d> &Lambda, T c) {
	IndexSet<Index1d> indexset;
	return indexset;
}

template <typename T>
IndexSet<Index1d>
C_periodic(const IndexSet<Index1d> &Lambda, T c) {
	IndexSet<Index1d> indexset;
	return indexset;
}


} // namespace lawa


#endif /* INDEXSET_TCC_ */
