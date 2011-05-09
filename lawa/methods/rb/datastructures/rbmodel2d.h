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
#include <lawa/methods/rb/solvers/solvers.h>
#include <lawa/operators/operator2d.h>
#include <lawa/righthandsides/rhs2d.h>
#include <lawa/settings/enum.h>

namespace lawa {

/* RBModel 2D: 
 *    This class contains all N-dependent data and functions needed for a reduced basis
 *    approximation.
 *    It contains a pointer to an associated Truth Model (which does not have to be initialized
 *    for online calculations).
 *
 *    The basis functions (\calN-dependent) are also stored here, so that we can reconstruct full 
 *    solutions without a complete truth model. 
 */

template <typename T, typename TruthModel>
class RBModel2D {

        typedef T (*theta_fctptr)(std::vector<T>& params); // Argumente -> eher auch RBThetaData-Objekt?
        typedef flens::GeMatrix<flens::FullStorage<T, cxxblas::ColMajor> >  FullColMatrixT;
        typedef flens::DenseVector<flens::Array<T> >                        DenseVectorT;  
        typedef Coefficients<Lexicographical,T,Index2D>                     CoeffVector;

    public:

    /* Public member functions */
        RBModel2D();

        void 
        set_current_param(const std::vector<T>& _param);

        std::vector<T>&
        get_current_param();

        unsigned int
        Q_a();

        unsigned int
        Q_f();

        unsigned int
        n_bf();

        void 
        attach_inner_product_op(Operator2D<T>& _inner_product_op);

        void
        set_truthmodel(TruthModel& _truthmodel);
                
    /* Public members */

        std::vector<theta_fctptr> theta_a;
        std::vector<theta_fctptr> theta_f;

        std::vector<FullColMatrixT>     RB_A_matrices;
        std::vector<DenseVectorT>       RB_F_vectors;
        std::vector<DenseVectorT>       RB_output_vectors;

        TruthModel*                     truth;

        std::vector<CoeffVector>        rb_basis_functions;

        void
        add_to_basis(const CoeffVector& sol);

        DenseVectorT
        RB_solve(unsigned int N, SolverCall call = call_cg);

        CoeffVector
        reconstruct_u_N(DenseVectorT u, unsigned int N);

    protected:

    /* Protected member functions */

        // Computes the inner product a(v,u) w.r.t. the operator a = inner_product_op
        T
        inner_product(const CoeffVector& v1, const CoeffVector& v2);

        // Update the (dense) RB matrices / vectors with the last basis function
        void
        update_RB_A_matrices();

        void
        update_RB_F_vectors();
        
        void
        update_RB_inner_product();
        
        void
        calculate_representors();

    /* Protected members */

        std::vector<T>  current_param;

        Operator2D<T>*  inner_product_op;
        
        FullColMatrixT  RB_inner_product;
        
        std::vector<FullColMatrixT> A_representors;
        std::vector<DenseVectorT>   F_representors;


};

} // namespace lawa


#include <lawa/methods/rb/datastructures/rbmodel2d.tcc>

#endif // LAWA_METHODS_RB_DATASTRUCTURES_RBMODEL2D_H
