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

#ifndef LAWA_ADAPTIVE_INDEXSET_H
#define LAWA_ADAPTIVE_INDEXSET_H 1

#include <lawa/math/math.h>
#include <lawa/support.h>

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

IndexSet<Index1D>
extractSpaceIndices(const IndexSet<Index2D> &Lambda)
{
	typedef IndexSet<Index2D>::const_iterator set2d_const_it;
	IndexSet<Index1D> ret(Lambda.d,Lambda.d_);
	for (set2d_const_it lambda=Lambda.begin(); lambda!=Lambda.end(); ++lambda) {
		ret.insert((*lambda).index2);
	}
	return ret;
}


template <typename T, DomainType Domain, Construction Cons>
IndexSet<Index1D>
C(const IndexSet<Index1D> &Lambda, T c, const Basis<T,Primal,Domain,Cons> &basis) {
    IndexSet<Index1D> ret(Lambda.d,Lambda.d_), tmp(Lambda.d, Lambda.d_);
    typedef typename IndexSet<Index1D>::const_iterator const_it;

    for (const_it lambda=Lambda.begin(); lambda!=Lambda.end(); ++lambda) {
        tmp = C((*lambda),c,basis);
        for (const_it mu=tmp.begin(); mu!=tmp.end(); ++mu) {
            if (Lambda.count(*mu) == 0) ret.insert(*mu);
        }
    }
    return ret;
}

template <typename T, DomainType Domain, Construction Cons>
IndexSet<Index1D>
C(const Index1D &lambda, T c, const Basis<T,Primal,Domain,Cons> &basis) {
    IndexSet<Index1D> ret(basis.d,basis.d_);
    const BSpline<T,Primal,Domain,Cons> phi(basis.mra,0);
    const Wavelet<T,Primal,Domain,Cons> psi(basis,0);

    C(lambda,c,basis.mra,basis,ret);
    return ret;
}


// Security zone interval
template <typename T, Construction Cons>
void
C(const Index1D &lambda, T c, const MRA<T,Primal,Interval,Cons> &mra,
  const Basis<T,Primal,Interval,Cons> &basis, IndexSet<Index1D> &ret)
{
    using std::min;
    using std::max;
    //ret.insert(Index1D(lambda.j,lambda.k,lambda.xtype));

    int j=lambda.j, jP1, k=lambda.k;
    XType xtype = lambda.xtype;
    int kMin_mra = basis.mra.rangeI(j).firstIndex(), kMax_mra = basis.mra.rangeI(j).lastIndex();
    if (lambda.xtype==XBSpline) {
        jP1=j;
        ret.insert(Index1D(j,std::max(k-2,kMin_mra),xtype));
        ret.insert(Index1D(j,std::max(k-1,kMin_mra),xtype));
        ret.insert(Index1D(j,std::min(k+1,kMax_mra),xtype));
        ret.insert(Index1D(j,std::min(k+2,kMax_mra),xtype));
    } else {
        jP1=j+1;
    }
    Support<T> supp;
    if (xtype==XBSpline)  supp=  mra.phi.support(j,k);
    else                  supp=basis.psi.support(j,k);

    T zLambda=0.5*(supp.l2+supp.l1);
    Support<T> contractedSupp(c*supp.l1 + (1-c)*zLambda, c*supp.l2 + (1-c)*zLambda);

    int kMin = basis.rangeJ(jP1).firstIndex(), kMax = basis.rangeJ(jP1).lastIndex();
    int kStart = std::min(std::max(iceil(contractedSupp.l1 * pow2i<T>(jP1)), kMin), kMax);
    assert((overlap(contractedSupp, basis.psi.support(jP1,kStart))>0));
    while ((kStart-1 >= kMin) && (overlap(contractedSupp, basis.psi.support(jP1,std::max(kStart-1, kMin)))>0)) {
        --kStart;
    }
    int kEnd = std::max(std::min(ifloor(contractedSupp.l2 * pow2i<T>(jP1)), kMax), kMin);
    assert((overlap(contractedSupp, basis.psi.support(jP1,kEnd))>0));
    while ((kEnd+1 <= kMax) && (overlap(contractedSupp, basis.psi.support(jP1,std::min(kEnd+1, kMax)))>0)) {
        ++kEnd;
    }

    for (int k=kStart; k<=kEnd; ++k) {
        ret.insert(Index1D(jP1,k,XWavelet));
    }
}

// Security zone periodic
template <typename T>
void
C(const Index1D &lambda, T c, const MRA<T,Primal,Periodic,CDF> &mra,
  const Basis<T,Primal,Periodic,CDF> &basis, IndexSet<Index1D> &ret)
{
    int j=lambda.j, k=lambda.k;
    XType xtype=lambda.xtype;
    //ret.insert(Index1D(j,k,xtype));
    if (xtype==XBSpline) {
        ret.insert(Index1D(j,(k-1 >= mra.rangeI(j).firstIndex()) ? k-1 : mra.rangeI(j).lastIndex()
                                    + ((1 - (mra.rangeI(j).firstIndex() - k+1))%mra.cardI(j)),xtype));
        ret.insert(Index1D(j,k+1 <= mra.rangeI(j).lastIndex() ? k+1 : mra.rangeI(j).firstIndex()
                                    - ((1 - (k+1 - mra.rangeI(j).lastIndex()))%mra.cardI(j)),xtype));
        ret.insert(Index1D(j,k-2 >= mra.rangeI(j).firstIndex() ? k-2 : mra.rangeI(j).lastIndex()
                                    + ((1 - (mra.rangeI(j).firstIndex() - k+2))%mra.cardI(j)),xtype));
        ret.insert(Index1D(j,k+2 <= mra.rangeI(j).lastIndex() ? k+2 : mra.rangeI(j).firstIndex()
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
                ret.insert(Index1D(j,k,XWavelet));
            }
        }
    }
    else {
        Support<T> contractedSupp, supp = basis.psi.psiR.support(j,k);
        T center = 0.5*(supp.l1 + supp.l2);
        contractedSupp = c*supp + (1-c)*center;
        /*    no wavelet indices on the same level?!
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
                ret.insert(Index1D(j,k,XWavelet));
            }
        }
        */

        long int kMin = floor( pow2i<T>(j+1)*contractedSupp.l1 - basis.psi.psiR.support(0,0).l2);
        long int kMax = ceil(pow2i<T>(j+1)*contractedSupp.l2 - basis.psi.psiR.support(0,0).l1);

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
                ret.insert(Index1D(j+1,k,XWavelet));
            }
        }

    }
}

// Security zone realline
template <typename T>
void
C(const Index1D &lambda, T c, const MRA<T,Primal,R,CDF> &mra,
  const Basis<T,Primal,R,CDF> &basis, IndexSet<Index1D> &ret)
{
    int j=lambda.j, k=lambda.k;
    XType xtype=lambda.xtype;
    if (xtype==XBSpline) {
        ret.insert(Index1D(j,k-1,xtype));
        ret.insert(Index1D(j,k+1,xtype));
        ret.insert(Index1D(j,k-2,xtype));
        ret.insert(Index1D(j,k+2,xtype));

        Support<T> contractedSupp, supp = mra.phi.support(j,k);
        T center = 0.5*(supp.l1 + supp.l2);
        contractedSupp.l1 = c*supp.l1 + (1-c)*center;
        contractedSupp.l2 = c*supp.l2 + (1-c)*center;

        int kMin = floor( pow2i<T>(j)*contractedSupp.l1 - basis.psi.support(0,0).l2);
        int kMax =  ceil( pow2i<T>(j)*contractedSupp.l2 - basis.psi.support(0,0).l1);
        for (int k1=kMin; k1<=kMax; ++k1) {
            if (overlap(contractedSupp, basis.psi.support(j,k1))>0) ret.insert(Index1D(j,k1,XWavelet));
        }
    }
    else {
        Support<T> contractedSupp, supp = basis.psi.support(j,k);
        T center = 0.5*(supp.l1 + supp.l2);
        contractedSupp.l1 = c*supp.l1 + (1-c)*center;
        contractedSupp.l2 = c*supp.l2 + (1-c)*center;
        /*
        long int kMin = floor( pow2i<T>(j)*contractedSupp.l1 - basis.psi.support(0,0).l2);
        long int kMax = ceil(pow2i<T>(j)*contractedSupp.l2 - basis.psi.support(0,0).l1);
        for (long int k1=kMin; k1<=kMax; ++k1) {
            if (overlap(contractedSupp, basis.psi.support(j,k1))>0) ret.insert(Index1D(j,k1,XWavelet));
        }
        */
        long int kMin = floor( pow2i<T>(j+1)*contractedSupp.l1 - basis.psi.support(0,0).l2);
        long int kMax = ceil(pow2i<T>(j+1)*contractedSupp.l2 - basis.psi.support(0,0).l1);
        for (long int k1=kMin; k1<=kMax; ++k1) {
            if (overlap(contractedSupp, basis.psi.support(j+1,k1))>0) ret.insert(Index1D(j+1,k1,XWavelet));
        }
    }
}


//Security zone 2D
template <typename T, typename Basis2D>
IndexSet<Index2D>
C(const IndexSet<Index2D> &Lambda, T c, const Basis2D &basis)
{
    typedef typename IndexSet<Index2D>::const_iterator const_it_2d;
    typedef typename IndexSet<Index1D>::const_iterator const_it;

    int d = Lambda.d, d_ = Lambda.d_;

    IndexSet<Index2D>  ret(d,d_);

    //Security zone of Lambda should not contain indices which are already in Lambda
    for (const_it_2d lambda=Lambda.begin(); lambda!=Lambda.end(); ++lambda) {
        IndexSet<Index1D > C_index1(d,d_), C_index2(d,d_);
        C_index1 = C((*lambda).index1, c, basis.first);
        C_index2 = C((*lambda).index2, c, basis.second);

        for (const_it it_C_index1=C_index1.begin(); it_C_index1!=C_index1.end(); ++it_C_index1) {
            if (Lambda.count(Index2D((*it_C_index1), (*lambda).index2))==0) {
                ret.insert(Index2D((*it_C_index1), (*lambda).index2));
            }
        }
        for (const_it it_C_index2=C_index2.begin(); it_C_index2!=C_index2.end(); ++it_C_index2) {
            if (Lambda.count(Index2D((*lambda).index1, (*it_C_index2)))==0) {
                ret.insert(Index2D((*lambda).index1, (*it_C_index2)));
            }
        }
/*
        for (const_it it_C_index1=C_index1.begin(); it_C_index1!=C_index1.end(); ++it_C_index1) {
            for (const_it it_C_index2=C_index2.begin(); it_C_index2!=C_index2.end(); ++it_C_index2) {
                if (Lambda.count(Index2D((*it_C_index1), (*it_C_index2)))==0) {
                    ret.insert(Index2D((*it_C_index1), (*it_C_index2)));
                }
            }
        }
*/
    }
    return ret;
}

template <typename T, typename Basis2D>
IndexSet<Index2D>
C_t(const IndexSet<Index2D> &Lambda, T c, const Basis2D &basis)
{
    typedef typename IndexSet<Index2D>::const_iterator const_it_2d;
    typedef typename IndexSet<Index1D>::const_iterator const_it;

    int d = Lambda.d, d_ = Lambda.d_;

    IndexSet<Index2D>  ret(d,d_);

    //Security zone of Lambda should not contain indices which are already in Lambda
    for (const_it_2d lambda=Lambda.begin(); lambda!=Lambda.end(); ++lambda) {
        IndexSet<Index1D > C_index1(d,d_), C_index2(d,d_);
        C_index1 = C((*lambda).index1, c, basis.first);
        //C_index2 = C((*lambda).index2, c, basis.second);

        for (const_it it_C_index1=C_index1.begin(); it_C_index1!=C_index1.end(); ++it_C_index1) {
            if (Lambda.count(Index2D((*it_C_index1), (*lambda).index2))==0) {
                ret.insert(Index2D((*it_C_index1), (*lambda).index2));
            }
        }
    }
    return ret;
}


//Security zone 3D
template <typename T, typename Basis3D>
IndexSet<Index3D>
C(const IndexSet<Index3D> &Lambda, T c, const Basis3D &basis)
{
    typedef typename IndexSet<Index3D>::const_iterator const_it_3d;
    typedef typename IndexSet<Index1D>::const_iterator const_it;

    int d = Lambda.d, d_ = Lambda.d_;

    IndexSet<Index3D>  ret(d,d_);

    //Security zone of Lambda should not contain indices which are already in Lambda
    for (const_it_3d lambda=Lambda.begin(); lambda!=Lambda.end(); ++lambda) {
        IndexSet<Index1D > C_index1(d,d_), C_index2(d,d_), C_index3(d,d_);
        C_index1 = C((*lambda).index1, c, basis.first);
        C_index2 = C((*lambda).index2, c, basis.second);
        C_index3 = C((*lambda).index3, c, basis.third);

        for (const_it it_C_index1=C_index1.begin(); it_C_index1!=C_index1.end(); ++it_C_index1) {
            if (Lambda.count(Index3D((*it_C_index1), (*lambda).index2, (*lambda).index3) )==0) {
                ret.insert(Index3D((*it_C_index1),   (*lambda).index2,  (*lambda).index3) );
            }
        }
        for (const_it it_C_index2=C_index2.begin(); it_C_index2!=C_index2.end(); ++it_C_index2) {
            if (Lambda.count(Index3D((*lambda).index1, (*it_C_index2), (*lambda).index3) )==0) {
                ret.insert(Index3D((*lambda).index1,   (*it_C_index2), (*lambda).index3) );
            }
        }
        for (const_it it_C_index3=C_index3.begin(); it_C_index3!=C_index3.end(); ++it_C_index3) {
            if (Lambda.count(Index3D((*lambda).index1, (*lambda).index2, (*it_C_index3) ) )==0) {
                ret.insert(Index3D((*lambda).index1,   (*lambda).index2, (*it_C_index3) ) );
            }
        }
    }
    return ret;
}

/*
 * Realizations of lambdaTilde1d
 */

template <typename T>
IndexSet<Index1D>
lambdaTilde1d_PDE(const Index1D &lambda, const Basis<T,Primal,R,CDF> &basis, int s_tilde, int jmin, int jmax, bool update)
{
	BSpline<T,Primal,R,CDF> phi(basis.mra);
	Wavelet<T,Primal,R,CDF> psi(basis);
    int j = lambda.j, k = lambda.k;
    int d = psi.d, d_= psi.d_;
    IndexSet<Index1D> ret(d,d_);
    Support<T> support_refbspline = phi.support(0,0);
    Support<T> support_refwavelet = psi.support(0,0);

    if (!update) {

        if (lambda.xtype == XBSpline) {
            Support<T> supp = phi.support(j,k);
            //cout << "lambdaTilde_R: Calculating IndexSet_R for BSpline with " << lambda << " " << " " << phi_col.singularSupport(j,k) << endl;

            // Inserting all indices corresponding to Bsplines with intersecting support using local compactness
            BSpline<T,Primal,R> phi_row(d);
            int kMin =  floor(pow2i<T>(j)*supp.l1 - support_refbspline.l2)-1;
            int kMax =   ceil(pow2i<T>(j)*supp.l2 - support_refbspline.l1)+1;
            for (int k_row=kMin; k_row<=kMax; ++k_row) {
                if (overlap(supp, phi.support(j,k_row)) > 0) {
                    //std::cout << "lambdaTilde: BSpline (" << j << ", " << k_row << "): " << phi_row.support(j,k_row) << " " << supp  << std::endl;
                    ret.insert(Index1D(j,k_row,XBSpline));
                }
            }

            // Inserting all indices corresponding to Wavelets with intersecting support using
            // a) local compactness  b) matrix compression  c) vanishing moments
            for (int j_row=j; j_row<=std::min(j+s_tilde, jmax); ++j_row) {        // realization of matrix compression via level threshold
                T Pow2i_Mjrow = pow2i<T>(-j_row);
                if (j_row>=j+2) {
                    DenseVector<Array<T> > singsupp = phi.singularSupport(j,k);
                    //cout << "LambdaTilde: Singular support phi_col = " << singpts;
                    for (int i=singsupp.firstIndex(); i<=singsupp.lastIndex(); ++i) {
                        int kMin = floor(pow2i<T>(j_row)*singsupp(i) - support_refwavelet.l2)-1;
                        int kMax =  ceil(pow2i<T>(j_row)*singsupp(i) - support_refwavelet.l1)+1;

                        for (int k_row=kMin; k_row<=kMax; ++k_row) {
                            Support<T> supp_row(Pow2i_Mjrow*(support_refwavelet.l1+k_row),Pow2i_Mjrow*(support_refwavelet.l2+k_row));// = psi.support(j_row,k_row);
                            if (((overlap(supp_row, supp) > 0)) && (!(distance(singsupp,supp_row) >= 0 ))) {
                                //std::cout << "LambdaTilde: Wavelet (" << j_row << ", k_row = " << k_row << "): " << psi.support(j_row,k_row) << " " << singsupp  << std::endl;
                                ret.insert(Index1D(j_row,k_row,XWavelet));
                            }
                        }
                    }
                }
                else {
                	int kMin = floor(pow2i<T>(j_row)*supp.l1 - support_refwavelet.l2)-1;
                	int kMax =  ceil(pow2i<T>(j_row)*supp.l2 - support_refwavelet.l1)+1;

                    for (int k_row=kMin; k_row<=kMax; ++k_row) {
                        Support<T> supp_row(Pow2i_Mjrow*(support_refwavelet.l1+k_row),Pow2i_Mjrow*(support_refwavelet.l2+k_row));// = psi.support(j_row,k_row);
                        if (overlap(supp, supp_row) > 0)  {
                            //std::cout << "LambdaTilde: Wavelet (" << j_row << ", k_row = " << k_row << "): " << psi.support(j_row,k_row) << " " << singsupp << std::endl;
                            ret.insert(Index1D(j_row,k_row,XWavelet));
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
                        ret.insert(Index1D(jmin,k_row,XBSpline));
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
                            if ((overlap(supp, supp_row) > 0) && (!(distance(singsupp,supp_row) >= 0 ))){
                                //std::cout << "LambdaTilde: Wavelet (" << j_row << ", k_row = " << k_row << "): " << psi.support(j_row,k_row) << " " << singsupp << std::endl;
                                ret.insert(Index1D(j_row,k_row,XWavelet));
                            }
                        }
                    }
                }
                else {
                    int kMin = floor(pow2i<T>(j_row)*supp.l1 - support_refwavelet.l2)-1;
                    int kMax =  ceil(pow2i<T>(j_row)*supp.l2 - support_refwavelet.l1)+1;

                    for (int k_row=kMin; k_row<=kMax; ++k_row) {
                        Support<T> supp_row(Pow2i_Mjrow*(support_refwavelet.l1+k_row),Pow2i_Mjrow*(support_refwavelet.l2+k_row));// = psi.support(j_row,k_row);
                        if ((overlap(supp, supp_row) > 0) && (!(distance(psi.singularSupport(j_row,k_row),supp) >= 0 ))) {
                            //std::cout << "LambdaTilde: Wavelet (" << j_row << ", k_row = " << k_row << "): " << psi.support(j_row,k_row) << " " << singsupp << std::endl;
                            ret.insert(Index1D(j_row,k_row,XWavelet));
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
                if (jmax >= j+level_diff) {                                                // level difference has to be large enough for vanishing entries due to vanishing moments
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

template <typename T>
IndexSet<Index1D>
lambdaTilde1d_PDE(const Index1D &lambda, const Basis<T,Primal,Periodic,CDF> &basis, int s_tilde, int jmin, int jmax, bool update)
{
	BSpline<T,Primal,Periodic,CDF> phi_col(basis.mra), phi_row(basis.mra);
	Wavelet<T,Primal,Periodic,CDF> psi_col(basis), psi_row(basis);
    int j = lambda.j, k = lambda.k;
    int d = basis.d, d_= basis.d_;
    IndexSet<Index1D> ret(d,d_);
    Support<T> support_refbspline = phi_col.phiR.support(0,0);
    Support<T> support_refwavelet = psi_col.psiR.support(0,0);

    if (!update) {

        if (lambda.xtype == XBSpline) {
            Support<T> supp = phi_col.phiR.support(j,k);

            // Inserting all indices corresponding to Bsplines with intersecting support using local compactness
            int kMin =  floor(pow2i<T>(j)*supp.l1 - support_refbspline.l2)-1;
            int kMax =   ceil(pow2i<T>(j)*supp.l2 - support_refbspline.l1)+1;
            for (int k_row=kMin; k_row<=kMax; ++k_row) {
                int k_row_per = k_row;
                if(k_row_per < basis.mra.rangeI(j).firstIndex()){
                	k_row_per = basis.mra.rangeI(j).lastIndex() + ((1 - (basis.mra.rangeI(j).firstIndex() - k_row_per))%basis.mra.cardI(j));
                }
                if(k_row_per > basis.mra.rangeI(j).lastIndex()){
                	k_row_per = basis.mra.rangeI(j).firstIndex() - ((1 - (k_row_per - basis.mra.rangeI(j).lastIndex()))%basis.mra.cardI(j));
                }
                if (overlap(phi_row.support(j,k_row_per), phi_col.support(j,k)) > 0) {
                	ret.insert(Index1D(jmin,k_row_per,XBSpline));
               }
            }

            // Inserting all indices corresponding to Wavelets with intersecting support using
            // a) local compactness  b) matrix compression  c) vanishing moments
            for (int j_row=j; j_row<=std::min(j+s_tilde, jmax); ++j_row) {        // realization of matrix compression via level threshold
                if (j_row>=j+2) {
                    DenseVector<Array<T> > singsupp = phi_col.phiR.singularSupport(j,k);
                    for (int i=singsupp.firstIndex(); i<=singsupp.lastIndex(); ++i) {
                        int kMin = floor(pow2i<T>(j_row)*singsupp(i) - support_refwavelet.l2)-1;
                        int kMax =  ceil(pow2i<T>(j_row)*singsupp(i) - support_refwavelet.l1)+1;

                        for (int k_row=kMin; k_row<=kMax; ++k_row) {
                            int k_row_per = k_row;
                            if(k_row_per < basis.rangeJ(j_row).firstIndex()){
                            	k_row_per = basis.rangeJ(j_row).lastIndex() + ((1 - (basis.rangeJ(j_row).firstIndex() - k_row_per))%basis.cardJ(j_row));
                           	}
                           	if(k_row_per> basis.rangeJ(j_row).lastIndex()){
                           		k_row_per = basis.rangeJ(j_row).firstIndex() - ((1 - (k_row_per - basis.rangeJ(j_row).lastIndex()))%basis.cardJ(j_row));
                           	}
                           	Support<T> supp_row = psi_row.support(j_row,k_row_per);
                            if (((overlap(supp_row, phi_col.support(j,k)) > 0)) &&
                            	(!(distance(phi_col.singularSupport(j,k),supp_row) >= 0 ))) {
                            	ret.insert(Index1D(j_row,k_row_per,XWavelet));
                            }
                        }
                    }
                }
                else {
                	int kMin = floor(pow2i<T>(j_row)*supp.l1 - support_refwavelet.l2)-1;
                	int kMax =  ceil(pow2i<T>(j_row)*supp.l2 - support_refwavelet.l1)+1;

                    for (int k_row=kMin; k_row<=kMax; ++k_row) {
                        int k_row_per = k_row;
                        if(k_row_per < basis.rangeJ(j_row).firstIndex()){
                            k_row_per = basis.rangeJ(j_row).lastIndex() + ((1 - (basis.rangeJ(j_row).firstIndex() - k_row_per))%basis.cardJ(j_row));
                        }
                        if(k_row_per > basis.rangeJ(j_row).lastIndex()){
                            k_row_per = basis.rangeJ(j_row).firstIndex() - ((1 - (k_row_per - basis.rangeJ(j_row).lastIndex()))%basis.cardJ(j_row));
                        }
                        if (overlap(psi_row.support(j_row,k_row_per), phi_col.support(j,k)) > 0) {
                            ret.insert(Index1D(j_row,k_row_per,XWavelet));
                        }
                    }
                }
            }
        }

        else {
            Support<T> supp = psi_col.support(j,k);

            // Inserting all indices corresponding to Bsplines with intersecting support using
            // a) local compactness  b) matrix compression  c) vanishing moments
            if (fabs(j - jmin) <= s_tilde) {
                int kMin = floor(pow2i<T>(j)*supp.l1 - support_refbspline.l2);
                int kMax =  ceil(pow2i<T>(j)*supp.l2 - support_refbspline.l1);
                for (int k_row=kMin; k_row<=kMax; ++k_row) {
                	int k_row_per = k_row;
                	if(k_row_per < basis.mra.rangeI(jmin).firstIndex()){
                	    k_row_per = basis.mra.rangeI(jmin).lastIndex() + ((1 - (basis.mra.rangeI(jmin).firstIndex() - k_row_per))%basis.mra.cardI(jmin));
                	}
                	if(k_row_per > basis.mra.rangeI(jmin).lastIndex()){
                	    k_row_per = basis.mra.rangeI(jmin).firstIndex() - ((1 - (k_row_per - basis.mra.rangeI(jmin).lastIndex()))%basis.mra.cardI(jmin));
                	}

                	if (overlap(phi_row.support(jmin,k_row_per), psi_col.support(j,k)) > 0) {
                	    ret.insert(Index1D(jmin,k_row_per,XBSpline));
                	}
                }
            }

			// Inserting all indices corresponding to Wavelets with intersecting support using
			// a) local compactness  b) matrix compression  c) vanishing moments
			for (int j_row=std::max(j-s_tilde,jmin); j_row<=std::min(j+s_tilde,jmax); ++j_row) {
				if (j_row>=j+2) {
					DenseVector<Array<T> > singsupp = psi_col.singularSupport(j,k);
					for (int i=singsupp.firstIndex(); i<=singsupp.lastIndex(); ++i) {
						int kMin = floor(pow2i<T>(j_row)*singsupp(i) - support_refwavelet.l2)-1;
						int kMax =  ceil(pow2i<T>(j_row)*singsupp(i) - support_refwavelet.l1)+1;

						for (int k_row=kMin; k_row<=kMax; ++k_row) {
							int k_row_per = k_row;
							if(k_row_per < basis.rangeJ(j_row).firstIndex()){
								k_row_per = basis.rangeJ(j_row).lastIndex() + ((1 - (basis.rangeJ(j_row).firstIndex() - k_row_per))%basis.cardJ(j_row));
							}
							if(k_row_per> basis.rangeJ(j_row).lastIndex()){
								k_row_per = basis.rangeJ(j_row).firstIndex() - ((1 - (k_row_per - basis.rangeJ(j_row).lastIndex()))%basis.cardJ(j_row));
							}
							Support<T> supp_row = psi_row.support(j_row,k_row_per);
							if (((overlap(supp_row, psi_col.support(j,k)) > 0)) &&
								(!(distance(psi_col.singularSupport(j,k),supp_row) >= 0 )) ) {
									ret.insert(Index1D(j_row,k_row_per,XWavelet));
							}
						}
					}
				}
				else {
					int kMin = floor(pow2i<T>(j_row)*supp.l1 - support_refwavelet.l2)-1;
					int kMax =  ceil(pow2i<T>(j_row)*supp.l2 - support_refwavelet.l1)+1;

					for (int k_row=kMin; k_row<=kMax; ++k_row) {
						int k_row_per = k_row;
						if(k_row_per < basis.rangeJ(j_row).firstIndex()){
							k_row_per = basis.rangeJ(j_row).lastIndex() + ((1 - (basis.rangeJ(j_row).firstIndex() - k_row_per))%basis.cardJ(j_row));
						}
						if(k_row_per > basis.rangeJ(j_row).lastIndex()){
							k_row_per = basis.rangeJ(j_row).firstIndex() - ((1 - (k_row_per - basis.rangeJ(j_row).lastIndex()))%basis.cardJ(j_row));
						}
						if ((overlap(psi_row.support(j_row,k_row_per), psi_col.support(j,k)) > 0) &&
						    !(distance(psi_row.singularSupport(j_row,k_row_per),supp) >= 0 ) ) {
						   ret.insert(Index1D(j_row,k_row_per,XWavelet));
						}
					}
				}
			}
        }
    }
    return ret;
}

template <typename T, Construction Cons>
IndexSet<Index1D>
lambdaTilde1d_PDE(const Index1D &lambda, const Basis<T,Primal,Interval,Cons> &basis, int s_tilde, int jmin, int jmax, bool update) {
	    using std::min;
	    using std::max;

	    BSpline<T,Primal,Interval,Cons> phi_col(basis.mra), phi_row(basis.mra);
	    Wavelet<T,Primal,Interval,Cons> psi_col(basis), psi_row(basis);
	    int j = lambda.j, k = lambda.k;
	    int d = basis.d, d_= basis.d_;
	    IndexSet<Index1D> ret(d,d_);

	    if (lambda.xtype==XBSpline) {

	    	Support<T> supp_col = phi_col.support(j,k);

	    	//Adding B-Splines
	        int kMin = basis.mra.rangeI(jmin).firstIndex(), kMax = basis.mra.rangeI(jmin).lastIndex();
	        int kStart = min(max(iceil(supp_col.l1 * pow2i<T>(jmin)),kMin), kMax);
	        //assert((overlap(supp_col, phi_row.support(jmin,kStart))>0));
	        while ((kStart-1 >= kMin) && (overlap(supp_col, phi_row.support(jmin,max(kStart-1, kMin)))>0)) {
	            --kStart;
	        }
	        int kEnd = max(min(ifloor(supp_col.l2 * pow2i<T>(jmin)),kMax), kMin);
	        assert((overlap(supp_col, phi_col.support(jmin,kEnd))>0));
	        while ((kEnd+1 <= kMax) && (overlap(supp_col, phi_row.support(jmin,min(kEnd+1,kMax)))>0)) {
	            ++kEnd;
	        }

	        for (int k_row=kStart; k_row<=kEnd; ++k_row) {
	            ret.insert(Index1D(jmin,k_row,XBSpline));
	        }

			//Adding Wavelets
	        for (int j_row=jmin; j_row<=min(jmin+s_tilde, jmax); ++j_row) {

	            int kMin = basis.rangeJ(j_row).firstIndex(), kMax = basis.rangeJ(j_row).lastIndex();
	            int kStart = min(max(iceil(supp_col.l1 * pow2i<T>(j_row)), kMin), kMax);
	            //assert((overlap(supp_col, psi_row.support(j_row,kStart))>0));
	            while (kStart-1>=kMin && overlap(supp_col,psi_row.support(j_row,max(kStart-1,kMin)))>0) {
	                --kStart;
	            }
	            int kEnd = max(min(ifloor(supp_col.l2 * pow2i<T>(j_row)), kMax), kMin);
	            //assert((overlap(supp_col, psi_row.support(j_row,kEnd))>0));
	            while (kEnd+1<=kMax && overlap(supp_col,psi_row.support(j_row,min(kEnd+1,kMax)))>0) {
	                ++kEnd;
	            }
	            for (int k_row=kStart; k_row<=kEnd; k_row++) {
	            	Range<int> rangeL = basis.rangeJL(j_row);
	            	Range<int> rangeR = basis.rangeJR(j_row);
	            	if ( (k_row <= rangeL.lastIndex()) || (k_row >= rangeR.firstIndex()) ) {
	            		ret.insert(Index1D(j_row,k_row,XWavelet));
	            		continue;
	            	}
	            	if  (!(distance(phi_col.singularSupport(j,k),psi_row.support(j_row,k_row)) >= 0 )) {	//singsupp!
	            		ret.insert(Index1D(j_row,k_row,XWavelet));
	            	}
	            }
	        }
	    }
	    else {
	    	Support<T> supp_col = psi_col.support(j,k);

	    	//Adding B-Splines
	        if (fabs(j - jmin) <= s_tilde) {
	            int kMin = basis.mra.rangeI(jmin).firstIndex(), kMax = basis.mra.rangeI(jmin).lastIndex();
	            int kStart = min(max(iceil(supp_col.l1 * pow2i<T>(jmin)), kMin), kMax);
	            //assert((overlap(supp_col, phi_row.support(jmin,kStart))>0));
	            while (kStart-1>=kMin && overlap(supp_col,phi_row.support(jmin,max(kStart-1,kMin)))>0) {
	                --kStart;
	            }
	            int kEnd = max(min(ifloor(supp_col.l2 * pow2i<T>(jmin)), kMax), kMin);
	            //assert((overlap(supp_col, phi_row.support(jmin,kEnd))>0));
	            while (kEnd+1<=kMax && overlap(supp_col,phi_row.support(jmin,min(kEnd+1,kMax)))>0) {
	                ++kEnd;
	            }

	            for (int k_row=kStart; k_row<=kEnd; ++k_row) {
	            	if  (distance(psi_col.singularSupport(j,k),phi_row.support(jmin,k_row)) < 0 ) {		//singsupp!
	            		ret.insert(Index1D(jmin,k_row,XBSpline));
	            	}
	            	else {
	            		if ( (k <= basis.rangeJL(j).lastIndex() )  ||
	            			 (k >= basis.rangeJR(j).firstIndex() )     ) {
								ret.insert(Index1D(jmin,k_row,XBSpline));
	            		}
	            	}
	            }
			}

	        //Adding Wavelets
	        for (int j_row=max(j-s_tilde,jmin); j_row<=min(j+s_tilde,jmax); ++j_row) {

	        	int kMin = basis.rangeJ(j_row).firstIndex(), kMax = basis.rangeJ(j_row).lastIndex();
	            int kStart = min(max(iceil(supp_col.l1 * pow2i<T>(j_row)), kMin), kMax);
	            //assert((overlap(supp_col, psi_row.support(j_row,kStart))>0));
	            while (kStart-1>=kMin && overlap(supp_col,psi_row.support(j_row,max(kStart-1,kMin)))>0) {
	                --kStart;
	            }
	            int kEnd = max(min(ifloor(supp_col.l2 * pow2i<T>(j_row)), kMax), kMin);
	            //assert((overlap(supp_col, psi_row.support(j_row,kEnd))>0));
	            while (kEnd+1<=kMax && overlap(supp_col,psi_row.support(j_row,min(kEnd+1,kMax)))>0) {
	                ++kEnd;
	            }
	            for (int k_row=kStart; k_row<=kEnd; ++k_row) {
	            	if (distance(psi_col.singularSupport(j,k),psi_row.support(j_row,k_row)) < 0 ) {
	            		ret.insert(Index1D(j_row,k_row,XWavelet));
	            		continue;
	            	}
	            	else {
	            		if ( (k_row <= basis.rangeJL(j_row).lastIndex() )  ||
	            			 (k_row >= basis.rangeJR(j_row).firstIndex() )     ) {
								ret.insert(Index1D(j_row,k_row,XWavelet));
	            		}
	            		continue;
	            	}

	            	if (distance(psi_row.singularSupport(j_row,k_row),supp_col) < 0 ) {
	            		ret.insert(Index1D(j_row,k_row,XWavelet));
	            	}
	            	else {
	            		if ( (k <= basis.rangeJL(j).lastIndex() )  ||
	            		     (k >= basis.rangeJR(j).firstIndex() )     ) {
								ret.insert(Index1D(j_row,k_row,XWavelet));
	            		}
	            	}

					/*
	            	if  ( (!(distance(psi_col.singularSupport(j,k),psi_row.support(j_row,k_row)) >= 0 )) &&	//singsupp!
	            		  (!(distance(psi_row.singularSupport(j_row,k_row), supp_col) >= 0 )) ) {			//singsupp!
	            		ret.insert(Index1D(j_row,k_row,XWavelet));
	                }
	                */
	            }
	        }
	    }
	    return ret;
}

} // namespace lawa

#endif // LAWA_ADAPTIVE_INDEXSET_H
