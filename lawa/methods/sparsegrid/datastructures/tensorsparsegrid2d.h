/*
  LAWA - Library for Adaptive Wavelet Applications.
  Copyright (C) 2008,2009 Sebastian Kestler, Mario Rometsch, Kristina Steih, Alexander Stippler.

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


#ifndef LAWA_METHODS_SPARSEGRID_DATASTRUCTURES_TENSORSPARSEGRID2D_H
#define LAWA_METHODS_SPARSEGRID_DATASTRUCTURES_TENSORSPARSEGRID2D_H 1

#include <lawa/methods/adaptive/datastructures/datastructures.h>
#include <lawa/methods/uniform/algorithms/blockassembler1d.h>
#include <lawa/operators/operators.h>


namespace lawa {

template <typename T, typename Basis2D, typename S1_x, typename S1_y, typename S2_x, typename S2_y>
class TensorSparseGrid2D {

    public:
        typedef flens::GeMatrix<flens::FullStorage<T, cxxblas::ColMajor> >  DenseMatrixT;
        typedef flens::SparseGeMatrix<flens::CRS<T,flens::CRS_General> >    SparseMatrixT;
        typedef flens::DenseVector<flens::Array<T> >                        DenseVectorT;

        typedef std::map<std::pair<int,int>, int ,lt_int_vs_int >           LevelPairMap;

        TensorSparseGrid2D(const Basis2D &_basis, const S1_x &_s1_x, const S1_y &_s1_y,
                           const S2_x &_s2_x, const S2_y &_s2_y,int _J, T _gamma);

        int
        getDimension() const;

        IndexSet<Index2D>
        getIndexSet() const;

        void
        toCoefficients(const DenseVectorT &vec,
                       Coefficients<Lexicographical,T,Index2D> &sparsegridcoefficients);

        DenseVectorT
        operator*(const DenseVectorT &v) const;

    private:
        DenseMatrixT
        block_multiplication(int i, int j, const DenseMatrixT &Xj) const;

        void
        assembleMatrices();

        const Basis2D       &basis;
        const S1_x          &s1_x;
        const S1_y          &s1_y;
        const S2_x          &s2_x;
        const S2_y          &s2_y;
        int                 j0_x, j0_y;
        int                 J;
        T                   gamma;
        BlockAssembler1D<T,typename Basis2D::FirstBasisType> blockassembler1d;

        int dim;
        std::vector<int*>                                 sg_blocks;
        LevelPairMap                                      levelpair_map;
        std::vector<SparseMatrixT>                        matrixblocks_s1_x, matrixblocks_s1_y,
                                                          matrixblocks_s2_x, matrixblocks_s2_y;
};

}   //namespace lawa

#include <lawa/methods/sparsegrid/datastructures/tensorsparsegrid2d.tcc>

#endif  // LAWA_METHODS_SPARSEGRID_DATASTRUCTURES_TENSORSPARSEGRID2D_H
