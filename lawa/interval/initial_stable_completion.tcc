/*
  LAWA - Library for Adaptive Wavelet Applications.
  Copyright (C) 2008,2009  Mario Rometsch, Alexander Stippler.

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

#include <lawa/aux/densify.h>
#include <extensions/flens/lapack_flens.h>

namespace lawa {

template <typename T, Construction ConsPrimal, Construction ConsDual>
void
initial_stable_completion(const MRA<T,Primal,Interval,ConsPrimal> &mra,
                          const MRA<T,Dual,Interval,ConsDual> &mra_,
                          GeMatrix<FullStorage<T,cxxblas::ColMajor> > &M1,
                          GeMatrix<FullStorage<T,cxxblas::ColMajor> > &M1_)
{
    typedef GeMatrix<FullStorage<T,cxxblas::ColMajor> > FullColMatrix;

    const RefinementMatrix<T,Interval,ConsPrimal> &M0 = mra.M0;
    M0.setLevel(mra_.min_j0);
    int d = mra.d, d_ = mra_.d_, l1 = mra.l1, l2 = mra.l2;

    FullColMatrix H(M0.numRows(), M0.numRows(), 
                    M0.firstRow(), M0.firstRow());
    FullColMatrix MM, MMTmp;
    densify(M0, MM, M0.firstRow(), M0.firstCol());
    FullColMatrix HHInv(M0.numRows(), M0.numRows(), 
                        M0.firstRow(), M0.firstRow()),
                  HHInvTmp;

    H.diag(0) = 1;
    HHInv.diag(0) = 1;
    for (int i=0; i<ifloor(mra.d/2.); ++i) {
        for (int p=MM.firstRow()+d-1-i; p<=MM.lastRow()-1-i; ++p) {
            if ((p-l1-i)%2==0) {
                H(p,p+1) = -MM(p,(p-l1-i)/2) / MM(p+1,(p-l1-i)/2);
            }
        }
        for (int p=MM.firstRow()+1+i; p<=MM.lastRow()-d+1+i; ++p) {
            if ((p-l2+i)%2==0) {
                H(p,p-1) = -MM(p,(p-l2+i)/2) / MM(p-1,(p-l2+i)/2);
            }
        }

        DenseVector<Array<int> > p(H.numRows());
        FullColMatrix HInv = H;
        trf(HInv,p);
        tri(HInv,p);

        blas::mm(NoTrans, NoTrans, 1., HHInv, HInv, 0., HHInvTmp);
        HHInv = HHInvTmp;

        blas::mm(NoTrans, NoTrans, 1., H, MM, 0., MMTmp);
        MM = MMTmp;
        H.diag(1) = 0;
        H.diag(-1) = 0;
    }

    FullColMatrix F(M0.numRows(), pow2i<T>(mra_.min_j0), M0.firstRow(), 0);
    if (d&1) {
        for (int q=F.firstCol(); q<=F.lastCol(); ++q) {
            if (q<pow2i<T>(mra_.min_j0-1)) {
                F(2*q,q) = MM(2*q+1,q);
                F(2*q+1,q) = -MM(2*q,q);
            } else {
                F(2*q,q) = -MM(2*q+1,q);
                F(2*q+1,q) = MM(2*q,q); 
            }
        }
    } else {
        for (int q=F.firstCol(); q<=F.lastCol(); ++q) {
            F(2*q+1,q) = 1.;
        }
    }

    //--- finally setting up M1 and M1_ ----------------------------------------
    //--- M1
    FullColMatrix Mj1initial;
    blas::mm(NoTrans,NoTrans, 1., HHInv, F, 0., Mj1initial);

    FullColMatrix Mj0, Mj0_;
    const RefinementMatrix<T,Interval,ConsDual> &M0_ = mra_.M0_;
    mra_.M0_.setLevel(mra_.min_j0);
    
    densify(M0, Mj0, M0.firstRow(), M0.firstCol());
    densify(M0_, Mj0_, M0.firstRow(), M0.firstCol());
    FullColMatrix Tmp(M0.numRows(), M0_.numRows());
    Tmp.diag(0) = 1;
    blas::mm(NoTrans,Trans, -.5, Mj0, Mj0_, 1., Tmp);
    blas::mm(NoTrans, NoTrans, 1., Tmp, Mj1initial, 0., M1);

    //--- M1_
    FullColMatrix Mj(pow2i<T>(mra_.min_j0+1)+l2-l1-1, 
                     pow2i<T>(mra_.min_j0+1)+l2-l1-1);
    Mj( _ , _(1,M0.numCols())) = Mj0;
    Mj( _ , _(M0.numCols()+1, Mj.lastCol())) = M1;

    DenseVector<Array<int> > p(Mj.numRows());
    FullColMatrix MjInv = Mj;
    trf(MjInv,p);
    tri(MjInv,p);

    FullColMatrix Mj1_Tmp = MjInv(_(pow2i<T>(mra_.min_j0)+d-1+1,
                                    pow2i<T>(mra_.min_j0+1)+d-1), _ );
    copy(Trans, Mj1_Tmp, M1_);
}

} // namespace lawa
