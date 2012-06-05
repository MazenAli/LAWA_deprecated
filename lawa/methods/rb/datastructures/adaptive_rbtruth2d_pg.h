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

#ifndef LAWA_METHODS_RB_DATASTRUCTURES_ADAPTIVE_RBTRUTH2D_PG_H
#define LAWA_METHODS_RB_DATASTRUCTURES_ADAPTIVE_RBTRUTH2D_PG_H 1

#include <lawa/methods/adaptive/datastructures/index.h>
#include <lawa/methods/adaptive/datastructures/indexset.h>
#include <lawa/methods/rb/datastructures/rbmodel2d.h>
#include <lawa/operators/operator2d.h>
#include <lawa/righthandsides/rhs2d.h>

namespace lawa {

/* Adaptive RBTruth 2D:
 *
 *    This class provides data and functions for an adaptive truth system, i.e. \calN-dependent
 *     members and methods, using adaptive operators and righthandsides.
 *
 *    It contains a pointer to the associated RB Model, as well as
 *    a pointer to an actual solver that is used for snapshot calculations.
 */
 
template <typename, typename> class RBModel2D;
 
template <typename T, typename TrialBasis, typename TestBasis, typename TrialPrec, typename TestPrec, typename TruthSolver, typename Compression>
class AdaptiveRBTruth2D_PG{

        typedef T (*theta_fctptr)(const std::vector<T>& params); // Argumente -> eher auch RBThetaData-Objekt?
        typedef flens::GeMatrix<flens::FullStorage<T, cxxblas::ColMajor> >  FullColMatrixT;
        typedef flens::SparseGeMatrix<flens::CRS<T,flens::CRS_General> >    SparseMatrixT;
        typedef flens::DenseVector<flens::Array<T> >                        DenseVectorT;  
        typedef Coefficients<Lexicographical,T,Index2D>                     CoeffVector;

    public:

        typedef TrialBasis                                         TrialBasisType;
        typedef TestBasis                                          TestBasisType;
        typedef TrialPrec                                          TrialPrecType;
        typedef TestPrec                                           TestPrecType;
    
    /* Public member functions */
                          
        AdaptiveRBTruth2D_PG(TrialBasis& _trialbasis, TestBasis& _testbasis, TrialPrec& _trialprec, TestPrec& _testprec, 
                          bool _use_inner_product = false, bool _use_A_matrix = false,
                          bool _orthwithX = false, bool _useinnprodX = false);
        
        void 
        attach_A_q (theta_fctptr theta_a_q, Operator2D<T>& A_q);
        
        void
        attach_A_q(Operator2D<T>& A_q);
        
        void
        attach_F_q(theta_fctptr theta_f_q, AdaptiveRhs<T, Index2D>& F_q);
        
        void
        attach_F_q(AdaptiveRhs<T, Index2D>& F_q);
        
        void 
        attach_test_inner_product_ops(AdaptiveOperator2D<T>& _test_inner_product_u_u_op, AdaptiveOperator2D<T>& _test_inner_product_v_v_op);
        
        void
        set_truthsolver(TruthSolver& _truthsolver);
        
        virtual void
        set_rb_model(RBModel2D<T, AdaptiveRBTruth2D_PG<T, TrialBasis, TestBasis, TrialPrec, TestPrec, TruthSolver, Compression> >& _rb);
        
        RBModel2D<T, AdaptiveRBTruth2D_PG<T, TrialBasis, TestBasis, TrialPrec, TestPrec, TruthSolver, Compression> >&
        get_rb_model();
        
        
        Coefficients<Lexicographical,T,Index2D>
        truth_solve();
        
        T
        get_trial_prec(const Index2D& index);
        
        T
        get_test_prec(const Index2D& index);
        
        void
        undo_trial_prec(CoeffVector& u_trial);
        
        void
        undo_test_prec(CoeffVector& u_test);
        
        virtual void
        add_new_basis_function(const CoeffVector& sol);
        
        void
        update_representors();
        
        void
        update_representor_norms();
        
        void
        write_riesz_representors(const std::string& directory_name = "offline_data/representors", int repr_nr = -1);
        
        void
        assemble_matrix(IndexSet<Index2D>& indexset_col, IndexSet<Index2D>& indexset_row, Operator2D<T>& op, SparseMatrixT& matrix);

        void
        assemble_vector(IndexSet<Index2D>& indexset, Operator2D<T>& op, DenseVectorT& vec);

        void
        assemble_all(IndexSet<Index2D>& indexset_col, IndexSet<Index2D>& indexset_row);

        void
        clear_assembled();

        /*
        void
        assemble_inner_product_matrix(IndexSet<Index2D>& trial_indexset, IndexSet<Index2D>& test_indexset);

        void
        assemble_A_operator_matrices(IndexSet<Index2D>& trial_indexset, IndexSet<Index2D>& test_indexset);
        */
            
        /* Wird dieses Ÿberhaupt gebraucht?
        T
        trial_inner_product(const CoeffVector& v1, const CoeffVector& v2);
    	*/

        /*
         * Inner Product in test space for two functions u in trial space
         */
        T
        test_inner_product_u_u(const CoeffVector& u1, const CoeffVector& u2);
        
        /*
         * Inner Product in test space for two functions v in test space
         */
        T
        test_inner_product_v_v(const CoeffVector& v1, const CoeffVector& v2);

        virtual T
        inner_product(const CoeffVector& v1, const CoeffVector& v2);

        T
        uncached_residual_dual_norm(const DenseVectorT& u_RB, const std::vector<T>& mu);

        T
        uncached_residual_dual_norm(const DenseVectorT& u_RB, const std::vector<T>& mu,
                                    Coefficients<Lexicographical,T,Index2D>& res_repr);
            
    /* Public members classes */
        // Wrapper class for affine structure on left hand side
       class Operator_LHS {

           private:

               AdaptiveRBTruth2D_PG<T, TrialBasis, TestBasis, TrialPrec, TestPrec, TruthSolver, Compression>* thisTruth;

           public:
               Operator_LHS(AdaptiveRBTruth2D_PG<T, TrialBasis, TestBasis, TrialPrec, TestPrec, TruthSolver, Compression>* _truth)
                   : thisTruth(_truth), compression(thisTruth->trialbasis), qa(-1){}

               T
               operator()(const Index2D &row_index, const Index2D &col_index);

               Coefficients<Lexicographical,T,Index2D>
               mv(const IndexSet<Index2D> &LambdaRow,
                  const Coefficients<Lexicographical,T,Index2D> &x);

               void
               toFlensSparseMatrix(const IndexSet<Index2D> &LambdaRow,
                                   const IndexSet<Index2D> &LambdaCol, SparseMatrixT &A, T tol);

               void
                apply(const Coefficients<Lexicographical,T,Index2D> &v, T eps,
                      const IndexSet<Index2D> &Lambda, Coefficients<Lexicographical,T,Index2D> &ret,
                      cxxblas::Transpose trans=cxxblas::NoTrans);

               Compression compression;

               int qa;
       };

        // Wrapper class for affine structure on right hand side
       class Operator_RHS {
           public:
               Operator_RHS(AdaptiveRBTruth2D_PG<T, TrialBasis, TestBasis, TrialPrec, TestPrec, TruthSolver, Compression>* _truth) : thisTruth(_truth){}

               T
               operator()(const Index2D &lambda);

               Coefficients<Lexicographical,T,Index2D>
               operator()(const IndexSet<Index2D> &Lambda);

               Coefficients<Lexicographical,T,Index2D>
               operator()(T tol);

           private:
               AdaptiveRBTruth2D_PG<T, TrialBasis, TestBasis, TrialPrec, TestPrec, TruthSolver, Compression>* thisTruth;
       };

       // Wrapper class for affine structure on left hand side
       // in the calculation of the Riesz representors
      class Operator_LHS_Representor {

          private:

              AdaptiveRBTruth2D_PG<T, TrialBasis, TestBasis, TrialPrec, TestPrec, TruthSolver, Compression>* thisTruth;

          public:
              Operator_LHS_Representor(AdaptiveRBTruth2D_PG<T, TrialBasis, TestBasis, TrialPrec, TestPrec, TruthSolver, Compression>* _truth)
                  : thisTruth(_truth), compression(thisTruth->trialbasis){}

              T
              operator()(const Index2D &row_index, const Index2D &col_index);

              Coefficients<Lexicographical,T,Index2D>
              mv(const IndexSet<Index2D> &LambdaRow,
                 const Coefficients<Lexicographical,T,Index2D> &x);

              void
              toFlensSparseMatrix(const IndexSet<Index2D> &LambdaRow,
                                  const IndexSet<Index2D> &LambdaCol, SparseMatrixT &A, T tol);

              void
              apply(const Coefficients<Lexicographical,T,Index2D> &v, T eps,
                    const IndexSet<Index2D> &Lambda, Coefficients<Lexicographical,T,Index2D> &ret,
                    cxxblas::Transpose trans=cxxblas::NoTrans);

              Compression compression;
      };

      // Wrapper class for affine structure on right hand side
      // in the calculation of Riesz representors for the bilinear
      // forms a^(q)
     class Operator_RHS_BilFormRepresentor {
         public:
             Operator_RHS_BilFormRepresentor(AdaptiveRBTruth2D_PG<T, TrialBasis, TestBasis, TrialPrec, TestPrec, TruthSolver, Compression>* _truth) : thisTruth(_truth){}

             T
             operator()(const Index2D &lambda);

             Coefficients<Lexicographical,T,Index2D>
             operator()(const IndexSet<Index2D> &Lambda);

             Coefficients<Lexicographical,T,Index2D>
             operator()(T tol);

             void
             set_current_op(Operator2D<T>& op);

             void
             set_current_bf(Coefficients<Lexicographical,T,Index2D>& bf);

         private:
             AdaptiveRBTruth2D_PG<T, TrialBasis, TestBasis, TrialPrec, TestPrec, TruthSolver, Compression>* thisTruth;

             Operator2D<T>*                              current_op;
             Coefficients<Lexicographical,T,Index2D>*    current_bf;
     };

      // Wrapper class for affine structure on right hand side
      // in the calculation of Riesz representors for the functional
      // forms f^(q)
     class Operator_RHS_FunctionalRepresentor {
         public:
             Operator_RHS_FunctionalRepresentor(AdaptiveRBTruth2D_PG<T, TrialBasis, TestBasis, TrialPrec, TestPrec, TruthSolver, Compression>* _truth) : thisTruth(_truth){}

             T
             operator()(const Index2D &lambda);

             Coefficients<Lexicographical,T,Index2D>
             operator()(const IndexSet<Index2D> &Lambda);

             Coefficients<Lexicographical,T,Index2D>
             operator()(T tol);

             void
             set_current_op(AdaptiveRhs<T, Index2D>& op);

         private:
             AdaptiveRBTruth2D_PG<T, TrialBasis, TestBasis, TrialPrec, TestPrec, TruthSolver, Compression>* thisTruth;

             AdaptiveRhs<T, Index2D>*                  current_op;
     };

     // Wrapper class for residual on right hand side
     // in the calculation of Riesz representor
     class Operator_Residual_Representor {

       public:
         Operator_Residual_Representor(AdaptiveRBTruth2D_PG<T, TrialBasis, TestBasis, TrialPrec, TestPrec, TruthSolver, Compression>* _truth,
                                       const std::vector<T>& mu, const DenseVectorT& _u_N)
           : thisTruth(_truth), eval_mu(mu), u_N(_u_N){};

         T
         operator()(const Index2D & lambda);

         Coefficients<Lexicographical,T,Index2D>
         operator()(const IndexSet<Index2D> &Lambda);

         Coefficients<Lexicographical,T,Index2D>
         operator()(T tol);

       private:
         AdaptiveRBTruth2D_PG<T, TrialBasis, TestBasis, TrialPrec, TestPrec, TruthSolver, Compression>* thisTruth;

         const std::vector<T>& eval_mu;
         const DenseVectorT& u_N;
     };


    /* Public members */
        
        TrialBasis&                             trialbasis, basis;
        TestBasis&                              testbasis;
        
        std::vector<Operator2D<T>*>             A_operators;
        std::vector<AdaptiveRhs<T, Index2D>*>   F_operators;
        // AdaptiveOperator2D<T>*                  trial_inner_product_op;
        AdaptiveOperator2D<T>*                  test_inner_product_u_u_op;
        AdaptiveOperator2D<T>*                  test_inner_product_v_v_op;
            
        TruthSolver*                            solver;
        
                        
        Operator_LHS lhs_op;
        Operator_RHS rhs_op;
                 
        Operator_LHS_Representor            repr_lhs_op;
        Operator_RHS_BilFormRepresentor     repr_rhs_A_op;
        Operator_RHS_FunctionalRepresentor  repr_rhs_F_op;
        
        bool use_inner_product_matrix;
        bool use_A_operator_matrices;
        
        bool assembled_inner_product_matrix;
        bool assembled_A_operator_matrices;
        bool assembled_prec_vec;

        //SparseMatrixT   trial_inner_product_matrix;
        SparseMatrixT   test_inner_product_u_u_matrix;
        SparseMatrixT   test_inner_product_v_v_matrix;
        DenseVectorT    test_prec_vec;
        std::vector<SparseMatrixT> A_operator_matrices;
      
    protected:
                
        RBModel2D<T, AdaptiveRBTruth2D_PG<T, TrialBasis, TestBasis, TrialPrec, TestPrec, TruthSolver, Compression> >*     rb;
        
        std::vector<CoeffVector>                F_representors; // Dim: 1 x Q_f
        std::vector<std::vector<CoeffVector> >  A_representors; // Dim: n x Q_a
        
        TrialPrec& trial_prec;
        TestPrec&  test_prec;
        Coefficients<Lexicographical, T, Index2D> trial_prec_data;
        Coefficients<Lexicographical, T, Index2D> test_prec_data;
        
        lt<AbsoluteValue,T> lt_obj;
        
        bool orth_wrt_X;
        bool use_innprod_X;
};
    
} // namespace lawa


#include <lawa/methods/rb/datastructures/adaptive_rbtruth2d_pg.tcc>

#endif // LAWA_METHODS_RB_DATASTRUCTURES_ADAPTIVE_RBTRUTH2D_PG_H
