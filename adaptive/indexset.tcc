/*
 * indexset.tcc
 *
 *  Created on: 10.08.2010
 *      Author: sebastian
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
C(const IndexSet<Index1d> &lambda, T c, const Basis<T,Primal,R,CDF> &basis) {
	IndexSet<Index1d> indexset;
	return indexset;
}

template <typename T, Construction Cons>
IndexSet<Index1d>
C(const IndexSet<Index1d> &lambda, T c, const Basis<T,Primal,Interval,Cons> &basis) {
	std::cout << "Security zone for interval basis!" << std::endl;
	IndexSet<Index1d> indexset;
	return indexset;
}

template <typename T>
IndexSet<Index1d>
C(const IndexSet<Index1d> &lambda, T c, const Basis<T,Primal,Periodic,CDF> &basis) {
	std::cout << "Security zone for periodic basis!" << std::endl;
	IndexSet<Index1d> indexset;
	return indexset;
}


} // namespace lawa


#endif /* INDEXSET_TCC_ */
