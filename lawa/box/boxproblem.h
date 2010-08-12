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
 
// TODO: change FullColMatrixT to SparseMatrix
/* TODO: enable sparse tensor basis 
 *      -> change J_x to some J_x(jy) in loops 
 *      -> add function J_x(int jy) = J_x to existing tensorbasis     
 */

#ifndef LAWA_BOX_BOXPROBLEM_H
#define LAWA_BOX_BOXPROBLEM_H 1

#include <lawa/box/boxindex.h>
#include <extensions/extensions.h>

namespace lawa{    

template<typename T, typename Basis>
class BoxProblem
{
    private:
        Basis basis;
   
    public: 
        BoxProblem(Basis _basis);
    
        template <typename BilinearForm>
        flens::SparseGeMatrix<flens::CRS<T,flens::CRS_General> >
        getStiffnessMatrix(BilinearForm& a, int J_x, int J_y, T tol = 10e-15);
    
        template <typename RHSIntegral>
        flens::DenseVector<flens::Array<T> >
        getRHS(RHSIntegral& rhs, int J_x, int J_y);
        
        template <typename Preconditioner>
        flens::DiagonalMatrix<T>    
        getPreconditioner(Preconditioner& P, int J_x, int J_y);
};

} // namespace lawa

#include <lawa/box/boxproblem.tcc>

#endif // LAWA_BOX_BOXPROBLEM_H