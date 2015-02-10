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
#include <lawa/methods/rb/deprecated/datastructures/datastructures.h>
#include <lawa/methods/rb/deprecated/solvers/solvers.h>
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

        typedef T (*theta_fctptr)(const std::vector<T>& params); // Argumente -> eher auch RBThetaData-Objekt?
        typedef flens::GeMatrix<flens::FullStorage<T, cxxblas::ColMajor> >  FullColMatrixT;
        typedef flens::DenseVector<flens::Array<T> >                        DenseVectorT;  
        typedef Coefficients<Lexicographical,T,Index2D>                     CoeffVector;

    public:
        typedef TruthModel      TruthModelType;

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
        Q_output();

        unsigned int
        n_bf();
        
        void
        attach_theta_a_q(theta_fctptr theta_a_q);
        
        void
        attach_theta_f_q(theta_fctptr theta_f_q);

        void
        attach_theta_output_q(theta_fctptr theta_output_q);

        void
        set_truthmodel(TruthModel& _truthmodel);
        
        void
        add_to_basis(const CoeffVector& sol);

        DenseVectorT
        RB_solve(unsigned int N, std::vector<T> param, SolverCall call = call_cg);
        
        DenseVectorT
        RB_solve(unsigned int N, SolverCall call = call_cg);

        CoeffVector
        reconstruct_u_N(DenseVectorT u, unsigned int N);
        
        T
        residual_dual_norm(const DenseVectorT& u_RB, const std::vector<T>& mu);
        
        // Lower bound for coercivity constant, min-Theta approach
        virtual T
        alpha_LB(std::vector<T>& _param);
        
        virtual T
        RB_errorbound(const DenseVectorT& u_RB, std::vector<T>& _param);

        void
        set_min_param(const std::vector<T>& _param);
    
        void
        set_max_param(const std::vector<T>& _param);
            
        void
        set_ref_param(const std::vector<T>& _param);
        
        void
        train_Greedy(const std::vector<T>& init_param, T tol, int Nmax, const char* filename = "Training.txt",
                     SolverCall call = call_cg, bool write_during_training = false, const char* foldername="training_adaptive");
        
        void
        train_strong_Greedy(const std::vector<T>& init_param, T tol, int Nmax, const char* filename = "Training.txt",
                     SolverCall call = call_cg, bool write_during_training = false, const char* foldername="training_adaptive");

        void
        generate_uniform_trainingset(std::vector<int>& param_nbs_per_dim);
        
        void
        generate_logarithmic_trainingset(std::vector<int>& param_nbs_per_dim);
        
        void
        generate_loglin2d_trainingset(std::vector<int>& param_nbs_per_dim);
        
        void
        write_basis_functions(const std::string& directory_name = "offline_data/bf", int bf_nr = -1);
        
        void
        read_basis_functions(const std::string& directory_name = "offline_data/bf");
        
        void
        write_RB_data(const std::string& directory_name = "offline_data");
        
        void
        read_RB_data(const std::string& directory_name = "offline_data");
                
    /* Public members */

        std::vector<theta_fctptr> theta_a;
        std::vector<theta_fctptr> theta_f;
        std::vector<theta_fctptr> theta_output;

        std::vector<FullColMatrixT>     RB_A_matrices;
        std::vector<DenseVectorT>       RB_F_vectors;
        std::vector<DenseVectorT>       RB_output_vectors;
        FullColMatrixT  				RB_inner_product;


        TruthModel*                     truth;

        std::vector<CoeffVector>        rb_basis_functions;
        
        FullColMatrixT                               F_F_representor_norms; // Speicherbedarf kann verringert werden..
        std::vector<FullColMatrixT>                  A_F_representor_norms;
        std::vector<std::vector<FullColMatrixT> >    A_A_representor_norms; //.. Ausnutzen der Symmetrie (Matrix als Vektor)
        FullColMatrixT                               output_output_representor_norms; // Speicherbedarf kann verringert werden..
        
        std::vector<std::vector<T> >  Xi_train;
        
        // Update the (dense) RB matrices / vectors with the last basis function
        void
        update_RB_A_matrices();

        void
        update_RB_F_vectors();
        
        void
        update_RB_inner_product();
        
    protected:

    /* Protected member functions */
      /*
        // Update the (dense) RB matrices / vectors with the last basis function
        void
        update_RB_A_matrices();

        void
        update_RB_F_vectors();
        
        void
        update_RB_inner_product();
        */

    /* Protected members */

        std::vector<T>  current_param;
        
        std::vector<T>  min_param;
        std::vector<T>  max_param;

        // Reference Parameter defining the inner product norm
        // Needed for min-Theta approach
        std::vector<T>  ref_param;

};

} // namespace lawa


#include <lawa/methods/rb/deprecated/datastructures/rbmodel2d.tcc>

#endif // LAWA_METHODS_RB_DATASTRUCTURES_RBMODEL2D_H
