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

#ifndef LAWA_METHODS_RB_DATASTRUCTURES_RBMODEL2D_H
#define LAWA_METHODS_RB_DATASTRUCTURES_RBMODEL2D_H 1

#include <vector>
#include <lawa/methods/adaptive/datastructures/index.h>
#include <lawa/methods/adaptive/datastructures/indexset.h>
#include <lawa/operators/operator2d.h>
#include <lawa/righthandsides/rhs2d.h>

namespace lawa {

/* RBModel 2D
 *
 */
 
template <typename T>
class RBModel2D {

    	typedef T (*theta_fctptr)(std::vector<T>& params); // Argumente -> eher auch RBThetaData-Objekt?
		typedef flens::GeMatrix<flens::FullStorage<T, cxxblas::ColMajor> >  FullColMatrixT;
		typedef flens::DenseVector<flens::Array<T> >                        DenseVectorT;  
        
	public:
        std::vector<theta_fctptr> theta_a;
        std::vector<theta_fctptr> theta_f;
                
        std::vector<FullColMatrixT> 	RB_A_matrices;
        std::vector<DenseVectorT> 		RB_F_vectors;
        std::vector<DenseVectorT>		RB_output_vectors;
        
        void 
        set_current_param(const std::vector<T>& _param);
        
        std::vector<T>& 
        get_current_param();
    	
    protected: 
    
    	RBModel2D();
        
        std::vector<T> current_param;
};
    
} // namespace lawa


#include <lawa/methods/rb/datastructures/rbmodel2d.tcc>

#endif // LAWA_METHODS_RB_DATASTRUCTURES_RBMODEL2D_H
