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

#ifndef LAWA_METHODS_RB_DATASTRUCTURES_UNIFORM_RBMODEL2D_H
#define LAWA_METHODS_RB_DATASTRUCTURES_UNIFORM_RBMODEL2D_H 1

#include <lawa/methods/adaptive/datastructures/index.h>
#include <lawa/methods/adaptive/datastructures/indexset.h>
#include <lawa/methods/rb/datastructures/rbmodel2d.h>
#include <lawa/operators/uniformoperator2d.h>
#include <lawa/preconditioners/nopreconditioner.h>
#include <lawa/righthandsides/rhs2d.h>

namespace lawa {

/* RBModel 2D
 *
 */
 
template <typename T, typename TruthSolver>
class UniformRBModel2D : public RBModel2D<T, TruthSolver> {

    	typedef T (*theta_fctptr)(std::vector<T>& params); // Argumente -> eher auch RBThetaData-Objekt?
		typedef flens::GeMatrix<flens::FullStorage<T, cxxblas::ColMajor> >  FullColMatrixT;
		typedef flens::DenseVector<flens::Array<T> >                        DenseVectorT;  
	
    public:

		/* Public member functions */
        
		UniformRBModel2D();
        
        void
        attach_A_q(theta_fctptr theta_a_q, UniformOperator2D<T>& A_q);
        
        void
        attach_F_q(theta_fctptr theta_f_q, Rhs2D<T>& F_q);
        
        void
        set_truthsolver(TruthSolver& _truthsolver);

        
        /* Public members */
        
        std::vector<UniformOperator2D<T>*> 	A_operators;
        std::vector<Rhs2D<T>*>			F_operators;
            	    
    	class Operator_LHS {
        	public:
            	Operator_LHS(UniformRBModel2D<T, TruthSolver>* _model) : thisModel(_model){}
                
				T
                operator()(XType row_xtype_x, int j1_x, int k1_x,
                           XType row_xtype_y, int j1_y, int k1_y,
                           XType col_xtype_x, int j2_x, int k2_x,
                           XType col_xtpye_y, int j2_y, int k2_y) const;
            
            private:
            	UniformRBModel2D<T, TruthSolver>* thisModel;
        };
        
        class Operator_RHS {
        	public:
            	Operator_RHS(UniformRBModel2D<T, TruthSolver>* _model) : thisModel(_model){}
                
				T
                operator()(XType xtype_x, int j_x, int k_x,
                           XType xtype_y, int j_y, int k_y) const;
            
            private:
            	UniformRBModel2D<T, TruthSolver>* thisModel;        
        };
        
        Operator_LHS 	lhs_op;
        Operator_RHS 	rhs_op;
    	
};
    
} // namespace lawa


#include <lawa/methods/rb/datastructures/uniform_rbmodel2d.tcc>

#endif // LAWA_METHODS_RB_DATASTRUCTURES_UNIFORM_RBMODEL2D_H
