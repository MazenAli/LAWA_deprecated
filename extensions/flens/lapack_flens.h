/*
 *   Copyright (c) 2007, Michael Lehn
 *
 *   All rights reserved.
 *
 *   Redistribution and use in source and binary forms, with or without
 *   modification, are permitted provided that the following conditions
 *   are met:
 *
 *   1) Redistributions of source code must retain the above copyright
 *      notice, this list of conditions and the following disclaimer.
 *   2) Redistributions in binary form must reproduce the above copyright
 *      notice, this list of conditions and the following disclaimer in
 *      the documentation and/or other materials provided with the
 *      distribution.
 *   3) Neither the name of the FLENS development group nor the names of
 *      its contributors may be used to endorse or promote products derived
 *      from this software without specific prior written permission.
 *
 *   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 *   "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 *   LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 *   A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 *   OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 *   SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 *   LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 *   DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 *   THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 *   (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 *   OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

#ifndef FLENS_LAPACK_FLENS_H
#define FLENS_LAPACK_FLENS_H 1

#include <complex>
#include <flens/vectortypes/impl/densevector.h>
#include <flens/matrixtypes/general/impl/gematrix.h>

namespace flens {

//== getrf ---------------------------------------------------------------------
template <typename FS>
    int
    trf(GeMatrix<FS> &A, DenseVector<Array<int> > &P);

//== getri ---------------------------------------------------------------------
template <typename FS>
    int
    tri(GeMatrix<FS> &A, DenseVector<Array<int> > &P);

//-- potrf ---------------------------------------------------------------------
template <typename FS>
    int
    trf(SyMatrix<FS> &A);

//== getrs ---------------------------------------------------------------------
template <typename MA, typename MB>
    int
    trs(cxxblas::Transpose trans, const GeMatrix<MA> &A,
        const DenseVector<Array<int> > &P, GeMatrix<MB> &B);

template <typename MA, typename VB>
    int
    trs(cxxblas::Transpose trans, const GeMatrix<MA> &A,
        const DenseVector<Array<int> > &P,
        DenseVector<VB> &B);

//== gesv ----------------------------------------------------------------------
template <typename MA, typename MB>
    int
    sv(GeMatrix<MA> &A, DenseVector<Array<int> > &P, GeMatrix<MB> &B);

template <typename MA, typename VB>
    int
    sv(GeMatrix<MA> &A, DenseVector<Array<int> > &P, DenseVector<VB> &B);

//== trtrs ---------------------------------------------------------------------
template <typename MA, typename MB>
    int
    trs(cxxblas::Transpose trans, const TrMatrix<MA> &A, GeMatrix<MB> &B);

template <typename MA, typename VB>
    int
    trs(cxxblas::Transpose trans, const TrMatrix<MA> &A, DenseVector<VB> &B);

//== geqrf ---------------------------------------------------------------------
template <typename MA, typename VT>
    int
    qrf(GeMatrix<MA> &A, DenseVector<VT> &tau);

//== geqp3 ---------------------------------------------------------------------
// QR decomposition with column pivoting
template <typename MA, typename VT>
	int
	qp3(GeMatrix<MA> &A, DenseVector<Array<int> >& p, DenseVector<VT> &tau);

//== orgqr ---------------------------------------------------------------------
template <typename MA, typename VT>
    int
    orgqr(GeMatrix<MA> &A, const DenseVector<VT> &tau);

//== ormqr ---------------------------------------------------------------------
template <typename MA, typename VT, typename MC>
    int
    ormqr(cxxblas::Side side, cxxblas::Transpose trans,
          const GeMatrix<MA> &A, const DenseVector<VT> &tau,
          GeMatrix<MC> &C);

//== gels ----------------------------------------------------------------------
template <typename MA, typename MB>
    int
    ls(cxxblas::Transpose trans, GeMatrix<MA> &A, GeMatrix<MB> &B);

//== gelss ---------------------------------------------------------------------
template <typename MA, typename MB>
    int
    lss(GeMatrix<MA> &A, GeMatrix<MB> &B);

//-- gees ----------------------------------------------------------------------
template <typename T>
struct ESSelectInfo
{
    typedef int (*pfunc)(T *vr, T *vi);
};

template <typename T>
struct ESSelectInfo<complex<T> >
{
    typedef int (*pfunc)(complex <T> *v);
};

// non complex version
template <typename MA, typename WR, typename WI, typename MVS>
    int
    es(bool CalcEv, bool sort,
       typename ESSelectInfo<typename MA::ElementType>::pfunc select,
       GeMatrix<MA> &A, int &sdim,
       DenseVector<WR> &wr, DenseVector<WI> &wi,
       GeMatrix<MVS> &VS);

// complex version
template <typename MA, typename VW, typename MVS>
    int
    es(bool CalcEv, bool sort,
       typename ESSelectInfo<typename MA::ElementType>::pfunc select,
       GeMatrix<MA> &A, int &sdim,
       DenseVector<VW> &w,
       GeMatrix<MVS> &VS);

//== geev,real -----------------------------------------------------------------
template <typename MA, typename WR, typename WI, typename VL, typename VR>
    int
    ev(bool leftEV, bool rightEV,
       GeMatrix<MA> &A, DenseVector<WR> &wr, DenseVector<WI> &wi,
       GeMatrix<VL> &vl, GeMatrix<VR> &vr);
       
//== ggev,real -----------------------------------------------------------------
template <typename MA, typename WR, typename WI, typename VL, typename VR>
    int
    gv(bool leftEV, bool rightEV,
       GeMatrix<MA> &A, GeMatrix<MA> &B, DenseVector<WR> &wr, DenseVector<WI> &wi,
       DenseVector<WR> &beta,
       GeMatrix<VL> &vl, GeMatrix<VR> &vr);

//== geev,complex --------------------------------------------------------------
template <typename MA, typename W, typename VL, typename VR>
    int
    ev(bool leftEv, bool rightEv,
       GeMatrix<MA> &A, DenseVector<W> &w, GeMatrix<VL> &vl, GeMatrix<VR> &vr);

//== syev ----------------------------------------------------------------------
template <typename MA, typename VW>
   int
   ev(bool compEV, SyMatrix<MA> &A, DenseVector<VW> &w);

//== heev ----------------------------------------------------------------------
template <typename MA, typename VW>
    int
    ev(bool compEV, HeMatrix<MA> &A, DenseVector<VW> &w);

//== gesvd ---------------------------------------------------------------------
template <typename MA, typename VS, typename VU, typename VVT>
    int
    svd(SVectorsJob jobu, SVectorsJob jobvt, GeMatrix<MA> &A,
        DenseVector<VS> &S, GeMatrix<VU> &U, GeMatrix<VVT> &VT);

// calls: svd(All,All,A,s,U,V)
template <typename MA, typename VS, typename MU, typename MV>
    int
    svd(GeMatrix<MA> &A, DenseVector<VS> &s, GeMatrix<MU> &U, GeMatrix<MV> &VT);
/*
//-- gesdd ---------------------------------------------------------------------
template <typename MA, typename VS, typename VU, typename VVT>
    int
    sdd(SVectorsJob jobz, GeMatrix<MA> &A,
        DenseVector<VS> &S, GeMatrix<VU> &U, GeMatrix<VVT> &VT);
*/
//------------------------------------------------------------------------------

} // namespace flens

#include <extensions/flens/lapack_flens.tcc>

#endif // FLENS_LAPACK_FLENS_H

