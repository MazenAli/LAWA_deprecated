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
IndexSet<Index>::IndexSet(int _d, int _d_)
: d(_d), d_(_d_)
{
}

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
void
C(const Index1d &lambda, T c, const BSpline<T,Primal,R,CDF> &phi, const Wavelet<T,Primal,R,CDF> &psi, IndexSet<Index1d> &ret) {
	int j=lambda.j, k=lambda.k;
	XType xtype=lambda.xtype;
	if (xtype==XBSpline) {
		ret.insert(Index1d(j,k,xtype));
		ret.insert(Index1d(j,k-1,xtype));
		ret.insert(Index1d(j,k+1,xtype));
		ret.insert(Index1d(j,k-2,xtype));
		ret.insert(Index1d(j,k+2,xtype));

		Support<T> contractedSupp, supp = phi.support(j,k);
		T center = 0.5*(supp.l1 + supp.l2);
		contractedSupp.l1 = c*supp.l1 + (1-c)*center;
		contractedSupp.l2 = c*supp.l2 + (1-c)*center;

		int kMin = floor( pow2i<T>(j)*contractedSupp.l1 - psi.support(0,0).l2);
		int kMax =  ceil( pow2i<T>(j)*contractedSupp.l2 - psi.support(0,0).l1);
		for (int k1=kMin; k1<=kMax; ++k1) {
			if (overlap(contractedSupp, psi.support(j,k1))>0) ret.insert(Index1d(j,k1,XWavelet));
		}
	}
	else {
		Support<T> contractedSupp, supp = psi.support(j,k);
		T center = 0.5*(supp.l1 + supp.l2);
		contractedSupp.l1 = c*supp.l1 + (1-c)*center;
		contractedSupp.l2 = c*supp.l2 + (1-c)*center;
		long int kMin = floor( pow2i<T>(j)*contractedSupp.l1 - psi.support(0,0).l2);
		long int kMax = ceil(pow2i<T>(j)*contractedSupp.l2 - psi.support(0,0).l1);
		for (long int k1=kMin; k1<=kMax; ++k1) {
			if (overlap(contractedSupp, psi.support(j,k1))>0) ret.insert(Index1d(j,k1,XWavelet));
		}
		kMin = floor( pow2i<T>(j+1)*contractedSupp.l1 - psi.support(0,0).l2);
		kMax = ceil(pow2i<T>(j+1)*contractedSupp.l2 - psi.support(0,0).l1);
		for (long int k1=kMin; k1<=kMax; ++k1) {
			if (overlap(contractedSupp, psi.support(j+1,k1))>0) ret.insert(Index1d(j+1,k1,XWavelet));
		}
	}
}

template <typename T>
void
C(const Index1d &lambda, T c, const MRA<T,Primal,Periodic,CDF> &mra, const Basis<T,Primal,Periodic,CDF> &basis, IndexSet<Index1d> &ret) {
	int j=lambda.j, k=lambda.k;
	XType xtype=lambda.xtype; 
	if (xtype==XBSpline) {
        std::cout << mra.rangeI(j) << std::endl;
		ret.insert(Index1d(j,k,xtype));
		ret.insert(Index1d(j,(k-1 >= mra.rangeI(j).firstIndex()) ? k-1 : mra.rangeI(j).lastIndex() 
		                            + ((1 - (mra.rangeI(j).firstIndex() - k+1))%mra.cardI(j)),xtype));
		ret.insert(Index1d(j,k+1 <= mra.rangeI(j).lastIndex() ? k+1 : mra.rangeI(j).firstIndex() 
		                            - ((1 - (k+1 - mra.rangeI(j).lastIndex()))%mra.cardI(j)),xtype));     
		ret.insert(Index1d(j,k-2 >= mra.rangeI(j).firstIndex() ? k-2 : mra.rangeI(j).lastIndex() 
		                            + ((1 - (mra.rangeI(j).firstIndex() - k+2))%mra.cardI(j)),xtype));
        ret.insert(Index1d(j,k+2 <= mra.rangeI(j).lastIndex() ? k+2 : mra.rangeI(j).firstIndex() 
                                    - ((1 - (k+2 - mra.rangeI(j).lastIndex()))%mra.cardI(j)),xtype));
		Support<T> contractedSupp, supp = mra.phi.phiR.support(j,k);
		T center = 0.5*(supp.l1 + supp.l2);
		contractedSupp = c*supp + (1-c)*center;	

		int kMin = floor( pow2i<T>(j)*contractedSupp.l1 - basis.psi.psiR.support(0,0).l2);
		int kMax =  ceil( pow2i<T>(j)*contractedSupp.l2 - basis.psi.psiR.support(0,0).l1);

		for (int k1=kMin; k1<=kMax; ++k1) {
			if (overlap(contractedSupp, basis.psi.psiR.support(j,k1))>0){
                int k = k1;
                if(k < basis.rangeJ(j).firstIndex()){
                    k = basis.rangeJ(j).lastIndex() + ((1 - (basis.rangeJ(j).firstIndex() - k))%basis.cardJ(j));
                }
                if(k > basis.rangeJ(j).lastIndex()){
                    k = basis.rangeJ(j).firstIndex() - ((1 - (k - basis.rangeJ(j).lastIndex()))%basis.cardJ(j));
                }
			    ret.insert(Index1d(j,k,XWavelet));
			}
		}
	}
	else {
		Support<T> contractedSupp, supp = basis.psi.psiR.support(j,k);
		T center = 0.5*(supp.l1 + supp.l2);
		contractedSupp = c*supp + (1-c)*center;

		long int kMin = floor( pow2i<T>(j)*contractedSupp.l1 - basis.psi.psiR.support(0,0).l2);
		long int kMax = ceil(pow2i<T>(j)*contractedSupp.l2 - basis.psi.psiR.support(0,0).l1);
        
		for (long int k1=kMin; k1<=kMax; ++k1) {
			if (overlap(contractedSupp, basis.psi.psiR.support(j,k1))>0){
			    int k = k1;
                if(k < basis.rangeJ(j).firstIndex()){
                    k = basis.rangeJ(j).lastIndex() + ((1 - (basis.rangeJ(j).firstIndex() - k))%basis.cardJ(j));
                }
                if(k > basis.rangeJ(j).lastIndex()){
                    k = basis.rangeJ(j).firstIndex() - ((1 - (k - basis.rangeJ(j).lastIndex()))%basis.cardJ(j));
                }
			    ret.insert(Index1d(j,k,XWavelet));
			}
		}
		
		kMin = floor( pow2i<T>(j+1)*contractedSupp.l1 - basis.psi.psiR.support(0,0).l2);
		kMax = ceil(pow2i<T>(j+1)*contractedSupp.l2 - basis.psi.psiR.support(0,0).l1);
        
		for (long int k1=kMin; k1<=kMax; ++k1) {
			if (overlap(contractedSupp, basis.psi.psiR.support(j+1,k1))>0)
			{
			    int k = k1;
                if(k < basis.rangeJ(j+1).firstIndex()){
                    k = basis.rangeJ(j+1).lastIndex() + ((1 - (basis.rangeJ(j+1).firstIndex() - k))%basis.cardJ(j+1));
                }
                if(k > basis.rangeJ(j+1).lastIndex()){
                    k = basis.rangeJ(j+1).firstIndex() - ((1 - (k - basis.rangeJ(j+1).lastIndex()))%basis.cardJ(j+1));
                }
			    ret.insert(Index1d(j+1,k,XWavelet));
			}
		}

	}
}

template <typename T>
IndexSet<Index1d>
C_realline(const IndexSet<Index1d> &Lambda, T c) {
	IndexSet<Index1d> ret(Lambda.d,Lambda.d_);
	typedef typename IndexSet<Index1d>::const_iterator const_it;

	const BSpline<T,Primal,R,CDF> phi(Lambda.d,0);
    const Wavelet<T,Primal,R,CDF> psi(Lambda.d,Lambda.d_,0);

    for (const_it lambda=Lambda.begin(); lambda!=Lambda.end(); ++lambda) {
		C((*lambda),c,phi,psi,ret);
    }
    return ret;
}

template <typename T>
IndexSet<Index1d>
C_interval(const IndexSet<Index1d> &Lambda, T c) {
	IndexSet<Index1d> indexset(Lambda.d,Lambda.d_);
	return indexset;
}

template <typename T>
IndexSet<Index1d>
C_periodic(const IndexSet<Index1d> &Lambda, T c) {
	IndexSet<Index1d> indexset(Lambda.d,Lambda.d_);
	typedef typename IndexSet<Index1d>::const_iterator const_it;

	const MRA<T,Primal,Periodic,CDF> mra(Lambda.d, Lambda.d_);
    const Basis<T,Primal,Periodic,CDF> basis(Lambda.d, Lambda.d_);

    for (const_it lambda=Lambda.begin(); lambda!=Lambda.end(); ++lambda) {
		C((*lambda), c, mra, basis, indexset);
    }
	return indexset;
}


} // namespace lawa


#endif /* INDEXSET_TCC_ */
