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

#ifndef LAWA_METHODS_RB_DATASTRUCTURES_ADAPTIVE_RBMODEL2D_H
#define LAWA_METHODS_RB_DATASTRUCTURES_ADAPTIVE_RBMODEL2D_H 1

#include <lawa/methods/adaptive/datastructures/index.h>
#include <lawa/methods/adaptive/datastructures/indexset.h>
#include <lawa/methods/rb/datastructures/rbmodel2d.h>
#include <lawa/operators/operator2d.h>
#include <lawa/righthandsides/rhs2d.h>

namespace lawa {

/* RBModel 2D
 *
 */
 
template <typename T, typename Basis, typename TruthSolver>
class AdaptiveRBModel2D : public RBModel2D<T, TruthSolver> {

    	typedef T (*theta_fctptr)(std::vector<T>& params); // Argumente -> eher auch RBThetaData-Objekt?
		typedef flens::GeMatrix<flens::FullStorage<T, cxxblas::ColMajor> >  FullColMatrixT;
		typedef flens::DenseVector<flens::Array<T> >                        DenseVectorT;  
        
	public:

		/* Public member functions */
        
		AdaptiveRBModel2D(Basis& _basis);
        
        void
        attach_A_q(theta_fctptr theta_a_q, Operator2D<T>& A_q);
        
        void
        attach_F_q(theta_fctptr theta_f_q, AdaptiveRhs<T, Index2D>& F_q);
        
		void
        set_truthsolver(TruthSolver& _truthsolver);
        
        /* Public members */
        
        Basis&										basis;
        
        std::vector<Operator2D<T>*> 	A_operators;
        std::vector<AdaptiveRhs<T, Index2D>*>		F_operators;
            	
    	class Operator_LHS {
        	
            typedef CompressionPDE2D<T, Basis> Compression;
        	
            private:
            	AdaptiveRBModel2D<T, Basis, TruthSolver>* thisModel;
                
        	public:
            	Operator_LHS(AdaptiveRBModel2D<T, Basis, TruthSolver>* _model)
                	: thisModel(_model), compression(thisModel->basis){}
                
				T
                operator()(const Index2D &row_index, const Index2D &col_index);
                
                Compression compression;
        };

        class Operator_RHS {
        	public:
            	Operator_RHS(AdaptiveRBModel2D<T, Basis, TruthSolver>* _model) : thisModel(_model){}
                
				T
                operator()(const Index2D &lambda);

                Coefficients<Lexicographical,T,Index2D>
                operator()(const IndexSet<Index2D> &Lambda);

                Coefficients<Lexicographical,T,Index2D>
                operator()(T tol);
            
            private:
            	AdaptiveRBModel2D<T, Basis, TruthSolver>* thisModel;        
        };
    	
        Operator_LHS lhs_op;
        Operator_RHS rhs_op;

};
    
} // namespace lawa


#include <lawa/methods/rb/datastructures/adaptive_rbmodel2d.tcc>

#endif // LAWA_METHODS_RB_DATASTRUCTURES_ADAPTIVE_RBMODEL2D_H
