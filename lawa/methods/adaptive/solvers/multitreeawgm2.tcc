#include <lawa/methods/adaptive/algorithms/algorithms.h>

#include <ios>
#include <iostream>
#include <fstream>
#include <iomanip>

namespace lawa {

template <typename Index, typename Basis, typename LocalOperator,
		  typename RHS, typename Preconditioner>
MultiTreeAWGM2<Index,Basis,LocalOperator,RHS,Preconditioner>::
MultiTreeAWGM2(const Basis &_basis, LocalOperator &_Op, RHS &_F, Preconditioner &_Prec)
: basis(_basis), Op(_Op), F(_F), Prec(_Prec), exact_sol(nullptr)
{}

template <typename Index, typename Basis, typename LocalOperator,
		  typename RHS, typename Preconditioner>
MultiTreeAWGM2<Index,Basis,LocalOperator,RHS,Preconditioner>::
MultiTreeAWGM2(const Basis &_basis, LocalOperator &_Op, RHS &_F, Preconditioner &_Prec,
				AWGM_Parameters& _awgm_params, IS_Parameters& _is_params)
: awgm_params(_awgm_params), is_params(_is_params),
  basis(_basis), Op(_Op), F(_F), Prec(_Prec), exact_sol(nullptr)
{}

template <typename Index, typename Basis, typename LocalOperator,
		  typename RHS, typename Preconditioner>
void
MultiTreeAWGM2<Index,Basis,LocalOperator,RHS,Preconditioner>::
solve(Coefficients<Lexicographical,T,Index> &u)
{
    IndexSet<Index> Lambda;
	if(u.size() > 0){
		Lambda = supp(u);
	}
	else{
		if(flens::IsSame<Index,Index2D>::value){
			getSparseGridIndexSet(basis,Lambda,2,0);
		}
		else{
			std::cerr << "Need some initial set in trial space! " << std::endl;
			exit(1);
		}
	}

	solve(u, Lambda);
}

template <typename Index, typename Basis, typename LocalOperator,
		  typename RHS, typename Preconditioner>
void
MultiTreeAWGM2<Index,Basis,LocalOperator,RHS,Preconditioner>::
solve(Coefficients<Lexicographical,T,Index> &u, IndexSet<Index>& Lambda)
{
    //---------------------------------------//
    //------- AWGM Initialization -----------//
    //---------------------------------------//

    /// Initialization of hash map vectors
    Coefficients<Lexicographical,T,Index2D> r(awgm_params.hashmapsize),
                                            p(awgm_params.hashmapsize),
                                            Ap(awgm_params.hashmapsize),
                                            res(awgm_params.hashmapsize),     // approximate residual for f-Au
                                            u_leafs(awgm_params.hashmapsize); // "leafs" of u

	FillWithZeros(Lambda,u_leafs);
	FillWithZeros(Lambda,res);

    // Default for initial cg tolerance computation if adaptive
    T res_norm = 1.;

    //---------------------------------------//
    //------- AWGM Iterations ---------------//
    //---------------------------------------//
    for(size_t awgm_its = 0; awgm_its <= awgm_params.max_its; ++awgm_its){

    	if (Lambda.size()>awgm_params.max_basissize){
    		break;
    	}

        awgm_info.sizeLambda.push_back(Lambda.size());

        //---------------------------------------//
        //------- CGLS  -------------------------//
        //---------------------------------------//

		// Re-initialize vectors from last iteration
        r.clear();
        p.clear();
        Ap.clear();
		FillWithZeros(Lambda,r);
		FillWithZeros(Lambda,p);
		FillWithZeros(Lambda,Ap);

		// CG Parameters
		T cg_tol;
		if(is_params.adaptive_tol){
			cg_tol = std::min(is_params.init_tol, is_params.res_reduction*res_norm);
		}
		else{
			cg_tol = is_params.absolute_tol;
		}

        if(awgm_params.verbose){
            std::cout << "******** Iteration " << awgm_its << std::endl;
            std::cout << std::right;
            std::cout << "   Current size of LambdaT " << std::setw(8) <<  Lambda.size() << std::endl;
            std::cout.precision();
            std::cout << std::left;
            std::cout << "   --- Starting CG with tolerance " << cg_tol << " ---" << std::endl;
        }

		// Local variables
		T alpha, beta, res_cg, res_cg_prev;
		T dummy=0.;

		// Initial step
		Op.eval(u,r,Prec);
		//r -= f;
	    for(auto& lambda : Lambda){
	    	r[lambda] -= Prec(lambda)*F(lambda);
	    }
		r *= -1;
		p = r;
		res_cg_prev = r*r;

		// CG Iterations
		for(size_t cg_its=0; cg_its <= is_params.max_its; ++cg_its){

			Ap.setToZero();						// Ap = A*p
			Op.eval(p,Ap,Prec);

			alpha = res_cg_prev / (p*Ap);
			u += alpha * p;
			r -= alpha * Ap;

			res_cg = r*r;

			if(is_params.verbose){
				std::cout.precision(12);
	            std::cout << "       CG Iteration " << std::setw(3) << cg_its << ": current error " << std::setw(18) << sqrt(res_cg) << std::endl;
			}

			if(std::sqrt(res_cg) <= cg_tol){
                std::cerr << "       CGLS stopped with error " << sqrt(res_cg) << " after "
                		  << cg_its << " iterations "<< std::endl;
                awgm_info.cg_its.push_back(cg_its);
				break;
			}

			beta  = res_cg/res_cg_prev;
			p *= beta;
			p += r;
			res_cg_prev = res_cg;
		}

		if(awgm_info.cg_its.size() < awgm_info.sizeLambda.size()){
            awgm_info.cg_its.push_back(is_params.max_its);
		}

        //---------------------------------------//
        //---- COMPUTE APPROXIMATE RESIDUAL -----//
        //---------------------------------------//

		// Compute cone around solution u
		IndexSet<Index2D>	Cdiff_u_leafs;	// new indizes in cone
        res.setToZero();
        extendMultiTree(basis, u_leafs, res, Cdiff_u_leafs, "standard", false, true);

        // Compute rhs on expanded test index set
        IndexSet<Index2D> LambdaRes = supp(res);

		// Compute residual Au - f on expanded test index set
        res.setToZero();
        Op.eval(u,res,Prec); 	// res = A*u
        //res -= f;
	    for(auto& lambda : LambdaRes){
	    	res[lambda] -= Prec(lambda)*F(lambda);
	    }
        res_norm = res.norm(2.);

        awgm_info.awgm_res.push_back(res_norm);
        awgm_info.sizeLambdaRes.push_back(LambdaRes.size());

        if(awgm_params.verbose){
            std::cout << "   --- Approximate Residual ---" << std::endl;
            std::cout << "       U_leafs:        " << std::setw(20) << u_leafs.size() << std::endl;
            std::cout << "       Residual:  " << std::setw(20) << res_norm
            								  << " (on " << std::setw(8) << std::right <<  LambdaRes.size()
            								  << std::left << " indizes)" << std::endl;
        }

        // Test for convergence
        if(res_norm <= awgm_params.tol){
            std::cerr << "AWGM target tolerance reached after " << awgm_its << " iterations: "
            		  << "Residual = " << res_norm << " "
                      << ", awgm_tol = " << awgm_params.tol << std::endl;

            if(awgm_params.print_info){
            	awgm_info.print();
            }

            if(awgm_params.plot_solution && flens::IsSame<Index,Index2D>::value){
            	if(awgm_params.verbose){
                    std::cout << "=====>  Plotting solution to file " << std::endl;
            	}
                plot2D<T,Basis,Preconditioner>(basis, u, Prec, exact_sol, 0., 1., 0., 1., 0.01, "awgm_cg_u_plot");
            }

            return;
        }

        //---------------------------------------//
        //---- COMPUTE NEXT INDEX SET -----------//
        //---------------------------------------//

        // Remove indizes from last iteration
        T P_Lambda_Res_square = 0.;
		for(auto& lambda : Lambda){
			P_Lambda_Res_square += std::pow(r[lambda], (T)2.);
			res.erase(lambda);
		}

		// Compute buckets
		T threshbound = std::sqrt(1.-awgm_params.alpha*awgm_params.alpha) * res.norm((T)2.)/std::sqrt(T(res.size()));
		Coefficients<Bucket,T,Index2D> r_bucket;
		r_bucket.bucketsort(res, threshbound);

        if(awgm_params.verbose){
            std::cout << "   --- Computing next index set ---" << std::endl;
            std::cout << "       Threshbound:         " << threshbound <<  std::endl;
            std::cout << "       ||P_{Lambda}r ||_2:  " << std::sqrt(P_Lambda_Res_square) << std::endl;
            std::cout << "       alpha*Residual:   " << awgm_params.alpha*res_norm << std::endl;
        }

		// Add buckets to dummy vector: p
		for (int i=0; i<(int)r_bucket.bucket_ell2norms.size(); ++i) {
			P_Lambda_Res_square += std::pow(r_bucket.bucket_ell2norms[i],2.0L);
			r_bucket.addBucketToCoefficients(p,i);

			if(awgm_params.verbose){
				std::cerr << "          Bucket " << i << ": L2-norm " << r_bucket.bucket_ell2norms[i] << std::endl;
			}

			if (P_Lambda_Res_square >= awgm_params.alpha*res_norm*awgm_params.alpha*res_norm) {
				break;
			}
		}

		// Replace indizes from last iteration
		for(auto& lambda : Lambda){
			res.insert(std::make_pair<decltype(lambda),T>(lambda,0.));
		}

		// New LambdaTrial: Add new indizes to u, complete to multitree
        IndexSet<Index2D> LambdaDiff;
        for(auto& coeff : p){
			if(u.find(coeff.first) == u.end()){
				completeMultiTree(basis, coeff.first, u, LambdaDiff, 0, true);
			}
		}

        // Compute new u_leafs
        u_leafs.clear();
        FillWithZeros(LambdaDiff,u_leafs);

        Lambda = supp(u);

        if(awgm_params.verbose){
            std::cout << "       Lambda:  raw extension  " <<  std::setw(8) << std::right <<  p.size() << std::endl;
            std::cout << "                multitree ext. " <<  std::setw(8) <<  LambdaDiff.size() << std::endl;
            std::cout << "                total size     " <<  std::setw(8) <<  Lambda.size() << std::endl;
            if(awgm_params.verbose_extra){
            	std::cout << LambdaDiff;
            }
            std::cout << std::endl;
        }
    }

    std::cerr << "AWGM reached maximal iteration number " << awgm_params.max_its << ": "
    		  << "Residual = " << res_norm << " "
              << ", awgm_tol = " << awgm_params.tol << std::endl;

    if(awgm_params.print_info){
    	if(awgm_params.verbose){
            std::cout << "=====>  Writing information to file " << std::endl;
    	}
    	awgm_info.print();
    }
    if(awgm_params.plot_solution && flens::IsSame<Index,Index2D>::value){
    	if(awgm_params.verbose){
            std::cout << "=====>  Plotting solution to file " << std::endl;
    	}
        plot2D<T,Basis,Preconditioner>(basis, u, Prec, exact_sol, 0., 1., 0., 1., 0.01, "awgm_cg_u_plot");
    }
}

template <typename Index, typename Basis, typename LocalOperator,
		  typename RHS, typename Preconditioner>
void
MultiTreeAWGM2<Index,Basis,LocalOperator,RHS,Preconditioner>::
set_sol(sol_fct_2d _sol)
{
	exact_sol = _sol;
}

} // namespace lawa