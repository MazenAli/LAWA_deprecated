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

#ifndef LAWA_METHODS_RB_DATASTRUCTURES_ADAPTIVE_RBTRUTH2D_H
#define LAWA_METHODS_RB_DATASTRUCTURES_ADAPTIVE_RBTRUTH2D_H 1

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
 
template <typename T, typename Basis, typename TruthSolver, typename Compression>
class AdaptiveRBTruth2D{

        typedef T (*theta_fctptr)(const std::vector<T>& params); // Argumente -> eher auch RBThetaData-Objekt?
        typedef flens::GeMatrix<flens::FullStorage<T, cxxblas::ColMajor> >  FullColMatrixT;
        typedef flens::SparseGeMatrix<flens::CRS<T,flens::CRS_General> >    SparseMatrixT;
        typedef flens::DenseVector<flens::Array<T> >                        DenseVectorT;  
        typedef Coefficients<Lexicographical,T,Index2D>                     CoeffVector;
        
        typedef CompressionWeightedPDE2D<T, Basis>                              WeightedCompression;
        typedef IndexsetTruthSolver<T, Basis, Index2D, WeightedCompression>     IndexsetSolver;

    public:

    /* Public member functions */
        
        AdaptiveRBTruth2D(Basis& _basis, bool _use_inner_product = false, bool _use_A_matrix = false,
                          bool _use_F_vector = false);
        
        void
        attach_A_q(theta_fctptr theta_a_q, Operator2D<T>& A_q);
        
        void
        attach_A_q(Operator2D<T>& A_q);
        
        void
        attach_F_q(theta_fctptr theta_f_q, AdaptiveRhs<T, Index2D>& F_q);
        
        void
        attach_F_q(AdaptiveRhs<T, Index2D>& F_q);
        
        void
        set_truthsolver(TruthSolver& _truthsolver);
        
        void
        set_rb_model(RBModel2D<T, AdaptiveRBTruth2D<T, Basis, TruthSolver, Compression> >& _rb);
        
        
        RBModel2D<T, AdaptiveRBTruth2D<T, Basis, TruthSolver, Compression> >&
        get_rb_model();
        
        void
        update_representors();
        
        void
        update_representor_norms();
        
        void
        write_riesz_representors(const std::string& directory_name = "offline_data/representors");
        
        void
        assemble_inner_product_matrix(IndexSet<Index2D>& indexset);

        void
        assemble_A_operator_matrices(IndexSet<Index2D>& indexset);
        
    /* Public members */
        
        Basis&                                  basis;
        
        std::vector<Operator2D<T>*>             A_operators;
        std::vector<AdaptiveRhs<T, Index2D>*>   F_operators;
        std::vector<AdaptiveRhs<T, Index2D>*>   output_operators;
        
        TruthSolver*                            solver;
                
         // Wrapper class for affine structure on left hand side       
        class Operator_LHS {
                        
            private:
                
                AdaptiveRBTruth2D<T, Basis, TruthSolver, Compression>* thisTruth;
                
            public:
                Operator_LHS(AdaptiveRBTruth2D<T, Basis, TruthSolver, Compression>* _truth)
                    : thisTruth(_truth), compression(thisTruth->basis), qa(-1){}
                
                T
                operator()(const Index2D &row_index, const Index2D &col_index);         
                
                Compression compression;
                
                int qa;
        };

         // Wrapper class for affine structure on right hand side       
        class Operator_RHS {
            public:
                Operator_RHS(AdaptiveRBTruth2D<T, Basis, TruthSolver, Compression>* _truth) : thisTruth(_truth){}
                
                T
                operator()(const Index2D &lambda);

                Coefficients<Lexicographical,T,Index2D>
                operator()(const IndexSet<Index2D> &Lambda);

                Coefficients<Lexicographical,T,Index2D>
                operator()(T tol);
            
            private:
                AdaptiveRBTruth2D<T, Basis, TruthSolver, Compression>* thisTruth;        
        };
        
        // Wrapper class for affine structure on output
       class Operator_output {
           public:
               Operator_output(AdaptiveRBTruth2D<T, Basis, TruthSolver, Compression>* _truth) : thisTruth(_truth){}

               T
               operator()(const Index2D &lambda);

               Coefficients<Lexicographical,T,Index2D>
               operator()(const IndexSet<Index2D> &Lambda);

               Coefficients<Lexicographical,T,Index2D>
               operator()(T tol);

           private:
               AdaptiveRBTruth2D<T, Basis, TruthSolver, Compression>* thisTruth;
       };

        Operator_LHS lhs_op;
        Operator_RHS rhs_op;
        Operator_output output_op;
         
          // Wrapper class for affine structure on left hand side
          // in the calculation of the Riesz representors       
         class Operator_LHS_Representor {

             private:

                 AdaptiveRBTruth2D<T, Basis, TruthSolver, Compression>* thisTruth;

             public:
                 Operator_LHS_Representor(AdaptiveRBTruth2D<T, Basis, TruthSolver, Compression>* _truth)
                     : thisTruth(_truth), compression(thisTruth->basis){}

                 T
                 operator()(const Index2D &row_index, const Index2D &col_index);

                 Compression compression;
         };
         
         // Wrapper class for affine structure on right hand side
         // in the calculation of Riesz representors for the bilinear
         // forms a^(q)               
        class Operator_RHS_BilFormRepresentor {
            public:
                Operator_RHS_BilFormRepresentor(AdaptiveRBTruth2D<T, Basis, TruthSolver, Compression>* _truth) : thisTruth(_truth){}
                
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
                AdaptiveRBTruth2D<T, Basis, TruthSolver, Compression>* thisTruth;
                
                Operator2D<T>*                              current_op;
                Coefficients<Lexicographical,T,Index2D>*    current_bf;
        };
        
         // Wrapper class for affine structure on right hand side
         // in the calculation of Riesz representors for the functional
         // forms f^(q)
        class Operator_RHS_FunctionalRepresentor {
            public:
                Operator_RHS_FunctionalRepresentor(AdaptiveRBTruth2D<T, Basis, TruthSolver, Compression>* _truth) : thisTruth(_truth){}
                
                T
                operator()(const Index2D &lambda);

                Coefficients<Lexicographical,T,Index2D>
                operator()(const IndexSet<Index2D> &Lambda);

                Coefficients<Lexicographical,T,Index2D>
                operator()(T tol);
                
                void
                set_current_op(AdaptiveRhs<T, Index2D>& op);
                
            private:
                AdaptiveRBTruth2D<T, Basis, TruthSolver, Compression>* thisTruth;
                
                AdaptiveRhs<T, Index2D>*                  current_op;
        };
        
        Operator_LHS_Representor            repr_lhs_op;
        Operator_RHS_BilFormRepresentor     repr_rhs_A_op;
        Operator_RHS_FunctionalRepresentor  repr_rhs_F_op;
        Operator_RHS_FunctionalRepresentor  repr_rhs_output_op;
        
        bool use_inner_product_matrix;
        bool use_A_operator_matrices;
        bool use_F_operator_vectors;
     //   bool use_output_operator_vectors;
        
        SparseMatrixT   inner_product_matrix;
        std::vector<SparseMatrixT> A_operator_matrices;
        
        
    private:
        
        RBModel2D<T, AdaptiveRBTruth2D<T, Basis, TruthSolver, Compression> >*     rb;
        
        std::vector<CoeffVector>                F_representors; // Dim: 1 x Q_f
        std::vector<std::vector<CoeffVector> >  A_representors; // Dim: n x Q_a
        std::vector<CoeffVector>                output_representors; // Dim: 1 x Q_m ??

};
    
} // namespace lawa


#include <lawa/methods/rb/datastructures/adaptive_rbtruth2d.tcc>

#endif // LAWA_METHODS_RB_DATASTRUCTURES_ADAPTIVE_RBTRUTH2D_H
