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
        	s << "  [" << (*lambda) << "], " << std::endl;
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

/*
 * Realizations of lambdaTilde1d
 */

template <typename T>
IndexSet<Index1d>
lambdaTilde1d_PDE(const Index1d &lambda, const BSpline<T,Primal,R,CDF> &phi, const Wavelet<T,Primal,R,CDF> &psi, int s_tilde, int jmin, int jmax, bool update)
{
	int j = lambda.j, k = lambda.k;
	int d = psi.d, d_= psi.d_;
	IndexSet<Index1d> ret(d,d_);
	Support<T> support_refbspline = phi.support(0,0);
	Support<T> support_refwavelet = psi.support(0,0);

	if (!update) {

		if (lambda.xtype == XBSpline) {
			Support<T> supp = phi.support(j,k);
			//cout << "lambdaTilde_R: Calculating IndexSet_R for BSpline with " << lambda << " " << " " << phi_col.singularSupport(j,k) << endl;

			// Inserting all indices corresponding to Bsplines with intersecting support using local compactness
			BSpline<T,Primal,R> phi_row(d);
			int kMin =  floor(pow2i<T>(j)*supp.l1 - support_refbspline.l2);
			int kMax =   ceil(pow2i<T>(j)*supp.l2 - support_refbspline.l1);
			for (int k_row=kMin; k_row<=kMax; ++k_row) {
				if (overlap(supp, phi.support(j,k_row)) > 0) {
					//std::cout << "lambdaTilde: BSpline (" << j << ", " << k_row << "): " << phi_row.support(j,k_row) << " " << supp  << std::endl;
					ret.insert(Index1d(j,k_row,XBSpline));
				}
			}

			// Inserting all indices corresponding to Wavelets with intersecting support using
			// a) local compactness  b) matrix compression  c) vanishing moments
			for (int j_row=j; j_row<=std::min(j+s_tilde, jmax); ++j_row) {		// realization of matrix compression via level threshold
				T Pow2i_Mjrow = pow2i<T>(-j_row);
				if (j_row>=j+2) {
					DenseVector<Array<T> > singsupp = phi.singularSupport(j,k);
					//cout << "LambdaTilde: Singular support phi_col = " << singpts;
					for (int i=singsupp.firstIndex(); i<=singsupp.lastIndex(); ++i) {
						int kMin = floor(pow2i<T>(j_row)*singsupp(i) - support_refwavelet.l2)-1;
						int kMax =  ceil(pow2i<T>(j_row)*singsupp(i) - support_refwavelet.l1)+1;

						for (int k_row=kMin; k_row<=kMax; ++k_row) {
							Support<T> supp_row(Pow2i_Mjrow*(support_refwavelet.l1+k_row),Pow2i_Mjrow*(support_refwavelet.l2+k_row));// = psi.support(j_row,k_row);
							if (((overlap(supp_row, supp) > 0)) && (!(distance(singsupp,supp_row) > 0 ))) {
								//std::cout << "LambdaTilde: Wavelet (" << j_row << ", k_row = " << k_row << "): " << psi.support(j_row,k_row) << " " << singsupp  << std::endl;
								ret.insert(Index1d(j_row,k_row,XWavelet));
							}
						}
					}
				}
				else {
					for (int k_row=kMin; k_row<=kMax; ++k_row) {
						Support<T> supp_row(Pow2i_Mjrow*(support_refwavelet.l1+k_row),Pow2i_Mjrow*(support_refwavelet.l2+k_row));// = psi.support(j_row,k_row);
						if (overlap(supp, supp_row) > 0)  {
							//std::cout << "LambdaTilde: Wavelet (" << j_row << ", k_row = " << k_row << "): " << psi.support(j_row,k_row) << " " << singsupp << std::endl;
							ret.insert(Index1d(j_row,k_row,XWavelet));
						}
					}
				}
			}
		}
		else {
			Support<T> supp = psi.support(j,k);
			//cout << "lambdaTilde_R: Calculating IndexSet_R for Wavelet with " << lambda << " " <<  psi_col.support(j,k) << " " << psi_col.singularSupport(j,k) << endl;

			// Inserting all indices corresponding to Bsplines with intersecting support using
			// a) local compactness  b) matrix compression  c) vanishing moments
			if (fabs(j - jmin) <= s_tilde) {
				int kMin = floor( pow2i<T>(jmin)*supp.l1 - phi.support(0,0).l2);
				int kMax =  ceil( pow2i<T>(jmin)*supp.l2 - phi.support(0,0).l1);
				for (int k_row=kMin; k_row<=kMax; ++k_row) {
					if (overlap(supp, phi.support(jmin,k_row)) > 0) {
						//std::cout << "lambdaTilde: BSpline (" << jmin << ", " << k_row << "): " << phi.support(jmin,k_row) << " " << supp  << std::endl;
						ret.insert(Index1d(jmin,k_row,XBSpline));
					}
				}
			}

			// Inserting all indices corresponding to Wavelets with intersecting support using
			// a) local compactness  b) matrix compression  c) vanishing moments
			for (int j_row=std::max(j-s_tilde,jmin); j_row<=std::min(j+s_tilde,jmax); ++j_row) {
				T Pow2i_Mjrow = pow2i<T>(-j_row);
				if (j_row>=j+2) {
					DenseVector<Array<T> > singsupp = psi.singularSupport(j,k);
					//cout << "LambdaTilde: Singular support psi_col_" << j << "," << k << " = " << singpts;
					for (int i=singsupp.firstIndex(); i<=singsupp.lastIndex(); ++i) {
						int kMin = floor(pow2i<T>(j_row)*singsupp(i) - support_refwavelet.l2)-1;
						int kMax =  ceil(pow2i<T>(j_row)*singsupp(i) - support_refwavelet.l1)+1;

						for (int k_row=kMin; k_row<=kMax; ++k_row) {
							Support<T> supp_row(Pow2i_Mjrow*(support_refwavelet.l1+k_row),Pow2i_Mjrow*(support_refwavelet.l2+k_row));// = psi.support(j_row,k_row);
							if ((overlap(supp, supp_row) > 0) && (!(distance(singsupp,supp_row) > 0 ))){
								//std::cout << "LambdaTilde: Wavelet (" << j_row << ", k_row = " << k_row << "): " << psi.support(j_row,k_row) << " " << singsupp << std::endl;
								ret.insert(Index1d(j_row,k_row,XWavelet));
							}
						}
					}
				}
				else {
					int kMin = floor(pow2i<T>(j_row)*supp.l1 - support_refwavelet.l2)-1;
					int kMax =  ceil(pow2i<T>(j_row)*supp.l2 - support_refwavelet.l1)+1;

					for (int k_row=kMin; k_row<=kMax; ++k_row) {
						Support<T> supp_row(Pow2i_Mjrow*(support_refwavelet.l1+k_row),Pow2i_Mjrow*(support_refwavelet.l2+k_row));// = psi.support(j_row,k_row);
						if (overlap(supp, supp_row) > 0) {
							//std::cout << "LambdaTilde: Wavelet (" << j_row << ", k_row = " << k_row << "): " << psi.support(j_row,k_row) << " " << singsupp << std::endl;
							ret.insert(Index1d(j_row,k_row,XWavelet));
						}
					}
				}
			}
		}
	}
/*
	// We assume that indices corresponding to non-zero values are already calculated up to level jmax-1 > jmin.
	// Therefore, we only have to add wavelets on level jmax.
	else {
		int j = lambda.j, k = lambda.k;
		int d = lambda.d, d_= lambda.d_;

		if (lambda.xtype == XBSpline) {
			assert(j == jmin);
			BSpline<T,Primal,R> phi_col(d);
			Support<T> supp = phi_col.support(j,k);
			//cout << "lambdaTilde_R: Calculating IndexSet_R for BSpline with " << lambda << " " <<  phi_col.support(j,k) << " " << phi_col.singularSupport(j,k) << endl;

			if (jmax <= j+s_tilde) {
				// Inserting all indices corresponding to Wavelets with intersecting support using
				// a) local compactness  b) matrix compression  c) vanishing moments
				Wavelet<T,Primal,R> psi_row(d,d_);
				DenseVector<Array<T> > singpts = phi_col.singularSupport(j,k);
				for (int i=singpts.firstIndex(); i<=singpts.lastIndex(); ++i) {
					int kMin =  ceil(pow2i(jmax)*singpts(i) - psi_row.support(0,0).l2);
					int kMax = floor(pow2i(jmax)*singpts(i) - psi_row.support(0,0).l1);
					for (int k_row=kMin; k_row<=kMax; ++k_row) {
						if ((overlap(psi_row.support(jmax,k_row), supp) == true)) {
							//cout << "LambdaTilde: kMin = " << kMin << ", kMax = " << kMax << ", k_row = " << k_row << ": " << psi_row.support(j_row,k_row) << endl;
							ret.insert(Index_R<T>(d,d_,jmax,k_row,XWavelet));
						}
					}
				}
			}
		}

		else {
			assert(j >= jmin);
			Wavelet<T,Primal,R> psi_col(d,d_);
			Support<T> supp = psi_col.support(j,k);
			//cout << "lambdaTilde_R: Calculating IndexSet_R for Wavelet with " << lambda << " " <<  psi_col.support(j,k) << " " << psi_col.singularSupport(j,k) << endl;

			// Inserting all indices corresponding to Wavelets with intersecting support using
			// a) local compactness  b) matrix compression  c) vanishing moments
			Wavelet<T,Primal,R> psi_row(d,d_);
			if (jmax <= j+s_tilde) {
				int level_diff = 3;
				if (d == 2) level_diff = 2;
				if (d == 3) level_diff = 4;
				if (jmax >= j+level_diff) {												// level difference has to be large enough for vanishing entries due to vanishing moments
					DenseVector<Array<T> > singpts = psi_col.singularSupport(j,k);
					//cout << "LambdaTilde: Singular support psi_col_" << j << "," << k << " = " << singpts;
					for (int i=singpts.firstIndex(); i<=singpts.lastIndex(); ++i) {
						int kMin =  ceil(pow2i(jmax)*singpts(i) - psi_row.support(0,0).l2);
						int kMax = floor(pow2i(jmax)*singpts(i) - psi_row.support(0,0).l1);
						for (int k_row=kMin; k_row<=kMax; ++k_row) {
							if ((overlap(psi_row.support(jmax,k_row), supp) == true)) {
								//cout << "LambdaTilde update: jmax = " << jmax << ", kMin = " << kMin << ", kMax = " << kMax << ", k_row = " << k_row << ": " << psi_row.support(jmax,k_row) << endl;
								ret.insert(Index_R<T>(d,d_,jmax,k_row,XWavelet));
							}
						}
					}
				}
				else {
					int kMin = ceil( pow2i(jmax)*supp.l1 - psi_row.support(0,0).l2);
					int kMax = floor(pow2i(jmax)*supp.l2 - psi_row.support(0,0).l1);
					for (int k_row=kMin; k_row<=kMax; ++k_row) {
						if ((overlap(supp, psi_row.support(jmax,k_row)) == true) && !((distSingularSupport(psi_row,psi_col,jmax,k_row,j,k) > 0.0 ))){
							//cout << "lambdaTilde_R: Wavelet (" << j_row << ", " << k_row << "): " << psi_row.support(j_row,k_row) << endl;
							ret.insert(Index_R<T>(d,d_,jmax,k_row,XWavelet));
						}
					}
				}
			}
		}


	}
	*/
	return ret;
}

} // namespace lawa


#endif /* INDEXSET_TCC_ */
