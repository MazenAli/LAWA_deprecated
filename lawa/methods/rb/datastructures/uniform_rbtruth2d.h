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

#ifndef LAWA_METHODS_RB_DATASTRUCTURES_UNIFORM_RBTRUTH2D_H
#define LAWA_METHODS_RB_DATASTRUCTURES_UNIFORM_RBTRUTH2D_H 1

#include <lawa/methods/adaptive/datastructures/index.h>
#include <lawa/methods/adaptive/datastructures/indexset.h>
#include <lawa/methods/rb/datastructures/rbmodel2d.h>
#include <lawa/operators/uniformoperator2d.h>
#include <lawa/preconditioners/nopreconditioner.h>
#include <lawa/righthandsides/rhs2d.h>
#include <lawa/settings/enum.h>

namespace lawa {

/* Uniform RBTruth 2D:
 *
 *    This class provides data and functions for an uniform truth system, i.e. \calN-dependent
 *     members and methods, using uniform operators and righthandsides.
 *
 *    It contains a pointer to the associated RB Model, as well as
 *    a pointer to an actual solver that is used for snapshot calculations.
 */

template <typename, typename> class RBModel2D;
 
template <typename T, typename TruthSolver>
class UniformRBTruth2D {

        typedef T (*theta_fctptr)(const std::vector<T>& params); // Argumente -> eher auch RBThetaData-Objekt?
        typedef flens::GeMatrix<flens::FullStorage<T, cxxblas::ColMajor> >  FullColMatrixT;
        typedef flens::DenseVector<flens::Array<T> >                        DenseVectorT;  
    
    public:

        static const VariationalFormulationType Type = Galerkin;
        typedef UniformOperator2D<T> OperatorType;

        /* Public member functions */
        
        UniformRBTruth2D();
        
        void
        attach_A_q(theta_fctptr theta_a_q, UniformOperator2D<T>& A_q);
        
        void
        attach_F_q(theta_fctptr theta_f_q, Rhs2D<T>& F_q);
        
        void
        set_truthsolver(TruthSolver& _truthsolver);
        
        void
        set_rb_model(RBModel2D<T, UniformRBTruth2D<T, TruthSolver> >& _rb);
        
        /* Public members */
        
        std::vector<UniformOperator2D<T>*>     A_operators;
        std::vector<Rhs2D<T>*>                F_operators;

        TruthSolver*                        solver;

         // Wrapper class for affine structure on left hand side       
        class Operator_LHS {
            public:
                Operator_LHS(UniformRBTruth2D<T, TruthSolver>* _truth) : thisTruth(_truth){}
                
                T
                operator()(XType row_xtype_x, int j1_x, int k1_x,
                           XType row_xtype_y, int j1_y, int k1_y,
                           XType col_xtype_x, int j2_x, int k2_x,
                           XType col_xtpye_y, int j2_y, int k2_y) const;
            
            private:
                UniformRBTruth2D<T, TruthSolver>* thisTruth;
        };
        
         // Wrapper class for affine structure on right hand side               
        class Operator_RHS {
            public:
                Operator_RHS(UniformRBTruth2D<T, TruthSolver>* _truth) : thisTruth(_truth){}
                
                T
                operator()(XType xtype_x, int j_x, int k_x,
                           XType xtype_y, int j_y, int k_y) const;
            
            private:
                UniformRBTruth2D<T, TruthSolver>* thisTruth;        
        };
        
        Operator_LHS     lhs_op;
        Operator_RHS     rhs_op;
        
    private:
    
        RBModel2D<T, UniformRBTruth2D<T, TruthSolver> >* rb;
        
};
    
} // namespace lawa


#include <lawa/methods/rb/datastructures/uniform_rbtruth2d.tcc>

#endif // LAWA_METHODS_RB_DATASTRUCTURES_UNIFORM_RBTRUTH2D_H
