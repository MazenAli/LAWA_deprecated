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

#ifndef EXTENSIONS_FLENS_SPARSE_BLAS_FLENS_H
#define EXTENSIONS_FLENS_SPARSE_BLAS_FLENS_H 1

#include <lawa/flensforlawa.h>
#include <extensions/flens/sparsematrix.h>
#include <extensions/flens/crs.h>

namespace flens {

// sparse_gemv
template <typename T, typename VX, typename VY>
void
mv(cxxblas::Transpose trans, T alpha, const SparseGeMatrix<flens::extensions::CRS<T> > &A,
   const DenseVector<VX> &x,
   typename DenseVector<VY>::ElementType beta, DenseVector<VY> &y);

// sparse_symv
template <typename T, CRS_Storage Storage, typename VX, typename VY>
void
mv(T alpha, const SparseSyMatrix<flens::extensions::CRS<T, Storage> > &A,
   const DenseVector<VX> &x,
   typename DenseVector<VY>::ElementType beta, DenseVector<VY> &y);

template <typename T>
void
my_mm_t(const SparseGeMatrix<flens::extensions::CRS<T> > &A,
        SparseGeMatrix<flens::extensions::CRS<T> > &B);

//   computes A^t A, requires A^t stored in CRS format
template <typename T>
void
my_mm_At_A(const SparseGeMatrix<flens::extensions::CRS<T> > &A,
           SparseGeMatrix<flens::extensions::CRS<T> > &B);

// sparse_gemm (sparse x dense)
template <typename T, CRS_Storage Storage>
void
mm(cxxblas::Transpose transA,
   cxxblas::Transpose transB,
   T alpha,
   const SparseGeMatrix<flens::extensions::CRS<T, Storage> > &A,
   const GeMatrix<FullStorage<T, cxxblas::ColMajor> > &B,
   T beta, GeMatrix<FullStorage<T, cxxblas::ColMajor> > &C);

} // flens

#include <extensions/flens/sparse_blas_flens.tcc>

#endif // EXTENSIONS_FLENS_SPARSE_BLAS_FLENS_H

