#include <lawa/methods/adaptive/algorithms/algorithms.h>

#include <ios>
#include <iostream>
#include <fstream>
#include <iomanip>

namespace lawa {

template <typename Index, typename TrialBasis, typename TestBasis,
		  typename LocalOperator, typename LocalOperatorTransp, typename RHS,
		  typename TrialPrec, typename TestPrec>
MultiTreeAWGM_PG<Index,TrialBasis,TestBasis,LocalOperator,LocalOperatorTransp,RHS,TrialPrec,TestPrec>::
MultiTreeAWGM_PG(const TrialBasis &_trialbasis, const TestBasis& _testbasis,
				LocalOperator &_Op, LocalOperatorTransp& _OpTransp, RHS &_F,
				TrialPrec &_trialPrec, TestPrec& _testPrec)
: trialbasis(_trialbasis), testbasis(_testbasis), Op(_Op), OpTransp(_OpTransp), F(_F),
  trialPrec(_trialPrec), testPrec(_testPrec), exact_sol(nullptr)
{}

template <typename Index, typename TrialBasis, typename TestBasis,
		  typename LocalOperator, typename LocalOperatorTransp, typename RHS,
		  typename TrialPrec, typename TestPrec>
MultiTreeAWGM_PG<Index,TrialBasis,TestBasis,LocalOperator,LocalOperatorTransp,RHS,TrialPrec,TestPrec>::
MultiTreeAWGM_PG(const TrialBasis &_trialbasis, const TestBasis& _testbasis,
				LocalOperator &_Op, LocalOperatorTransp& _OpTransp, RHS &_F,
				TrialPrec &_trialPrec, TestPrec& _testPrec,
				AWGM_PG_Parameters& _awgm_params, IS_Parameters& _is_params)
: awgm_params(_awgm_params), is_params(_is_params),
  trialbasis(_trialbasis), testbasis(_testbasis), Op(_Op), OpTransp(_OpTransp), F(_F),
  trialPrec(_trialPrec), testPrec(_testPrec), exact_sol(nullptr)
{}

template <typename Index, typename TrialBasis, typename TestBasis,
		  typename LocalOperator, typename LocalOperatorTransp, typename RHS,
		  typename TrialPrec, typename TestPrec>
RHS&
MultiTreeAWGM_PG<Index,TrialBasis,TestBasis,LocalOperator,LocalOperatorTransp,RHS,TrialPrec,TestPrec>::
get_rhs()
{
	return F;
}

template <typename Index, typename TrialBasis, typename TestBasis,
		  typename LocalOperator, typename LocalOperatorTransp, typename RHS,
		  typename TrialPrec, typename TestPrec>
LocalOperator&
MultiTreeAWGM_PG<Index,TrialBasis,TestBasis,LocalOperator,LocalOperatorTransp,RHS,TrialPrec,TestPrec>::
get_lhs()
{
	return Op;
}

template <typename Index, typename TrialBasis, typename TestBasis,
		  typename LocalOperator, typename LocalOperatorTransp, typename RHS,
		  typename TrialPrec, typename TestPrec>
const TrialBasis&
MultiTreeAWGM_PG<Index,TrialBasis,TestBasis,LocalOperator,LocalOperatorTransp,RHS,TrialPrec,TestPrec>::
get_trialbasis()
{
	return trialbasis;
}

template <typename Index, typename TrialBasis, typename TestBasis,
		  typename LocalOperator, typename LocalOperatorTransp, typename RHS,
		  typename TrialPrec, typename TestPrec>
const TestBasis&
MultiTreeAWGM_PG<Index,TrialBasis,TestBasis,LocalOperator,LocalOperatorTransp,RHS,TrialPrec,TestPrec>::
get_testbasis()
{
	return testbasis;
}

template <typename Index, typename TrialBasis, typename TestBasis,
		  typename LocalOperator, typename LocalOperatorTransp, typename RHS,
		  typename TrialPrec, typename TestPrec>
void
MultiTreeAWGM_PG<Index,TrialBasis,TestBasis,LocalOperator,LocalOperatorTransp,RHS,TrialPrec,TestPrec>::
set_initial_indexsets(const IndexSet<Index> _LambdaTrial, const IndexSet<Index> _LambdaTest)
{
	default_init_Lambda_trial = _LambdaTrial;
	default_init_Lambda_test  = _LambdaTest;
}

template <typename Index, typename TrialBasis, typename TestBasis,
		  typename LocalOperator, typename LocalOperatorTransp, typename RHS,
		  typename TrialPrec, typename TestPrec>
void
MultiTreeAWGM_PG<Index,TrialBasis,TestBasis,LocalOperator,LocalOperatorTransp,RHS,TrialPrec,TestPrec>::
remove_preconditioner(Coefficients<Lexicographical,T,Index> &u)
{
	for(auto& el : u){
		el.second *= trialPrec(el.first);
	}
}


template <typename Index, typename TrialBasis, typename TestBasis,
		  typename LocalOperator, typename LocalOperatorTransp, typename RHS,
		  typename TrialPrec, typename TestPrec>
Coefficients<Lexicographical,typename LocalOperator::T,Index>
MultiTreeAWGM_PG<Index,TrialBasis,TestBasis,LocalOperator,LocalOperatorTransp,RHS,TrialPrec,TestPrec>::
solve()
{
    IndexSet<Index> LambdaTrial = default_init_Lambda_trial;
    IndexSet<Index> LambdaTest = default_init_Lambda_test;

	Coefficients<Lexicographical,T,Index> u;
	FillWithZeros(LambdaTrial,u);

	solve(u,LambdaTrial,LambdaTest);

	return u;
}

template <typename Index, typename TrialBasis, typename TestBasis,
		  typename LocalOperator, typename LocalOperatorTransp, typename RHS,
		  typename TrialPrec, typename TestPrec>
void
MultiTreeAWGM_PG<Index,TrialBasis,TestBasis,LocalOperator,LocalOperatorTransp,RHS,TrialPrec,TestPrec>::
solve(Coefficients<Lexicographical,T,Index> &u)
{
    IndexSet<Index> LambdaTrial, LambdaTest;
	if(u.size() > 0){

		// Take support of u as initial index set
		LambdaTrial = supp(u);

		// Check if this is the default trial index set
		bool is_default_set = true;
		if(LambdaTrial.size() != default_init_Lambda_trial){
			is_default_set = false;
		}
		else{
			for(auto& lambda : default_init_Lambda_trial){
				if(LambdaTrial.find(lambda) == LambdaTrial.end()){
					is_default_set = false;
					break;
				}
			}
		}

		// If it is not the default set, we have to get some stable
		// test index set -> stable expansion
		if(!is_default_set){
			std::cerr << "Computing stable expansion as test set. " << std::endl;
			Coefficients<Lexicographical,T,Index2D> Lambda_aux;
			getStableExpansion(trialbasis, testbasis, LambdaTrial, Lambda_aux);
			LambdaTest = supp(Lambda_aux);
		}
		// Else we just use the default test set
		else{
			LambdaTest = default_init_Lambda_test;
		}
	}
	else{
		LambdaTrial = default_init_Lambda_trial;
		LambdaTest = default_init_Lambda_test;
		FillWithZeros(LambdaTrial,u);
	}

	solve(u, LambdaTrial, LambdaTest);
}

template <typename Index, typename TrialBasis, typename TestBasis,
		  typename LocalOperator, typename LocalOperatorTransp, typename RHS,
		  typename TrialPrec, typename TestPrec>
void
MultiTreeAWGM_PG<Index,TrialBasis,TestBasis,LocalOperator,LocalOperatorTransp,RHS,TrialPrec,TestPrec>::
solve(Coefficients<Lexicographical,T,Index> &u, IndexSet<Index>& LambdaTrial, IndexSet<Index>& LambdaTest)
{
    //---------------------------------------//
    //------- AWGM Initialization -----------//
    //---------------------------------------//

    /// Initialization of hash map vectors
    Coefficients<Lexicographical,T,Index2D> r(awgm_params.hashmapsize_test),
    										s(awgm_params.hashmapsize_trial),
                                            p(awgm_params.hashmapsize_trial),
                                            Ap(awgm_params.hashmapsize_test),
                                            res(awgm_params.hashmapsize_test),     // approximate residual for f-Au
                                    		resNE(awgm_params.hashmapsize_test),   // cone around u_leafs
                                            u_leafs(awgm_params.hashmapsize_trial); // "leafs" of u

	FillWithZeros(LambdaTrial,u_leafs);
	FillWithZeros(LambdaTest,res);

    // Default for initial cgls tolerance computation if adaptive
    T resNE_norm = 1.;
    T res_norm = 1.;

    //---------------------------------------//
    //------- AWGM Iterations ---------------//
    //---------------------------------------//
    for(size_t awgm_its = 0; awgm_its <= awgm_params.max_its; ++awgm_its){

    	if (LambdaTrial.size()>awgm_params.max_basissize){
    		break;
    	}

        awgm_info.sizeLambdaTrial.push_back(LambdaTrial.size());
        awgm_info.sizeLambdaTest.push_back(LambdaTest.size());

        //---------------------------------------//
        //------- CGLS  -------------------------//
        //---------------------------------------//

		// Re-initialize vectors from last iteration
        r.clear();
        s.clear();
        p.clear();
        Ap.clear();
		FillWithZeros(LambdaTest,r);
		FillWithZeros(LambdaTrial,p);
		FillWithZeros(LambdaTrial,s);
		FillWithZeros(LambdaTest,Ap);

		// CGLS Parameters
		T cgls_tol;
		if(is_params.adaptive_tol){
			cgls_tol = std::min(is_params.init_tol, is_params.res_reduction*resNE_norm);
		}
		else{
			cgls_tol = is_params.absolute_tol;
		}

        if(awgm_params.verbose){
            std::cout << "******** Iteration " << awgm_its << std::endl;
            std::cout << std::right;
            std::cout << "   Current size of LambdaTrial " << std::setw(8) <<  LambdaTrial.size() << std::endl;
            std::cout << "   Current size of LambdaTest  " << std::setw(8) << LambdaTest.size() << std::endl << std::endl;
            std::cout.precision();
            std::cout << std::left;
            std::cout << "   --- Starting CGLS with tolerance " << cgls_tol << " ---" << std::endl;
        }

		// Local variables
		T alpha, beta, gamma_cgls, gamma_cgls_Prev, res_cgls;
		T dummy=0.;

		// Initial step
		Op.eval(u,r,trialPrec,testPrec);
		//r -= f;
	    for(auto& lambda : LambdaTest){
	    	r[lambda] -= testPrec(lambda)*F(lambda);
	    }
		r *= -1;
		OpTransp.eval(r,s,testPrec,trialPrec);
		p = s;
		gamma_cgls_Prev = s*s;

		// CGLS Iterations
		for(size_t cgls_its=0; cgls_its <= is_params.max_its; ++cgls_its){

			Ap.setToZero();						// Ap = A*p
			Op.eval(p,Ap,trialPrec,testPrec);

			alpha = gamma_cgls_Prev / (Ap*Ap);
			u += alpha * p;
			r -= alpha * Ap;

			s.setToZero();
			OpTransp.eval(r,s,testPrec,trialPrec);	// s = A^T*r

			gamma_cgls = s*s;
			res_cgls = r.norm(2.);

			if(is_params.verbose){
				std::cout.precision(12);
	            std::cout << "       CGLS Iteration " << std::setw(3) << cgls_its << ": current error NE " << std::setw(18) << sqrt(gamma_cgls)
	            												 << ", Au-f " << std::setw(18) << res_cgls << std::endl;
			}

			if(std::sqrt(gamma_cgls) <= cgls_tol){
                std::cerr << "       CGLS stopped with NE error " << sqrt(gamma_cgls) << ", error Au=f = " << res_cgls
                	 << " after " << cgls_its << " iterations "<< std::endl;
                awgm_info.cgls_its.push_back(cgls_its);
				break;
			}

			beta  = gamma_cgls/gamma_cgls_Prev;
			p *= beta;
			p += s;
			gamma_cgls_Prev = gamma_cgls;
		}

		if(awgm_info.cgls_its.size() < awgm_info.sizeLambdaTrial.size()){
            awgm_info.cgls_its.push_back(is_params.max_its);
		}

        //---------------------------------------//
        //---- COMPUTE APPROXIMATE RESIDUAL -----//
        //---------------------------------------//

		// Compute cone around solution u
		IndexSet<Index2D>	Cdiff_u_leafs;	// new indizes in cone
        res.setToZero();

        extendMultiTree(trialbasis, u_leafs, resNE, Cdiff_u_leafs, "standard", false, true);

        // Compute stable expansion in test basis
        IndexSet<Index2D> LambdaResNE = supp(resNE);
        if(awgm_params.reset_res){
        	res.clear();
        }
        getStableExpansion(trialbasis, testbasis, LambdaResNE, res);

        // Compute rhs on expanded test index set
        IndexSet<Index2D> LambdaRes = supp(res);

		// Compute residual Au - f on expanded test index set
        res.setToZero();
        Op.eval(u,res,trialPrec,testPrec); 	// res = A*u
        //res -= f;
	    for(auto& lambda : LambdaRes){
	    	res[lambda] -= testPrec(lambda)*F(lambda);
	    }
        res_norm = res.norm(2.);

        awgm_info.awgm_res.push_back(res_norm);
        awgm_info.sizeLambdaRes.push_back(LambdaRes.size());

        // Compute residual of NE: A^TA*u - A^T*f
        resNE.setToZero();
        OpTransp.eval(res,resNE,testPrec,trialPrec); 	// resNE = A^T*res
        resNE_norm = resNE.norm(2.);

        awgm_info.awgm_resNE.push_back(resNE_norm);
        awgm_info.sizeLambdaResNE.push_back(LambdaResNE.size());

        if(awgm_params.verbose){
            std::cout << "   --- Approximate Residual ---" << std::endl;
            std::cout << "       U_leafs:        " << std::setw(20) << u_leafs.size() << std::endl;
            std::cout << "       Residual Au-f:  " << std::setw(20) << res_norm
            									   << " (on " << std::setw(8) << std::right <<  LambdaRes.size() << std::left << " indizes)" << std::endl;
            std::cout << "       Residual NE  :  " << std::setw(20) << resNE_norm
            									   << " (on " << std::setw(8) << std::right << LambdaResNE.size()<< std::left << " indizes)" << std::endl;
        }

        // Test for convergence
        if(resNE_norm <= awgm_params.tol){
            std::cerr << "AWGM target tolerance reached after " << awgm_its << " iterations: "
            		  << "Residual NE = " << resNE_norm << " "
            		  << ", Residual Au-f = " << res_norm << " "
                      << ", awgm_tol = " << awgm_params.tol << std::endl<< std::endl;

            if(awgm_params.print_info){
            	if(awgm_params.verbose){
                    std::cout << "=====>  Writing information to file " << std::endl << std::endl;
            	}
            	awgm_info.print(awgm_params.info_filename.c_str());
            }

            if(awgm_params.plot_solution && flens::IsSame<Index,Index2D>::value){
            	if(awgm_params.verbose){
                    std::cout << "=====>  Plotting solution to file " << std::endl << std::endl;
            	}
                plot2D<T,TrialBasis,TrialPrec>(trialbasis, u, trialPrec, exact_sol, 0., 1., 0., 1., 0.01, awgm_params.plot_filename.c_str());
            }

            return;
        }

        //---------------------------------------//
        //---- COMPUTE NEXT INDEX SET -----------//
        //---------------------------------------//

        // Remove indizes from last iteration
        T P_Lambda_ResNE_square = 0.;
		for(auto& lambda : LambdaTrial){
			P_Lambda_ResNE_square += std::pow(s[lambda], (T)2.);
			resNE.erase(lambda);
		}

		// Compute buckets
		T threshbound = std::sqrt(1.-awgm_params.alpha*awgm_params.alpha) * resNE.norm((T)2.)/std::sqrt(T(resNE.size()));
		Coefficients<Bucket,T,Index2D> r_bucket;
		r_bucket.bucketsort(resNE, threshbound);

        if(awgm_params.verbose){
            std::cout << "   --- Computing next index set ---" << std::endl;
            std::cout << "       Threshbound:         " << threshbound <<  std::endl;
            std::cout << "       ||P_{Lambda}r ||_2:  " << std::sqrt(P_Lambda_ResNE_square) << std::endl;
            std::cout << "       alpha*Residual_NE:   " << awgm_params.alpha*resNE_norm << std::endl;
        }

		// Add buckets to dummy vector: p
		for (int i=0; i<(int)r_bucket.bucket_ell2norms.size(); ++i) {
			P_Lambda_ResNE_square += std::pow(r_bucket.bucket_ell2norms[i],2.0L);
			r_bucket.addBucketToCoefficients(p,i);

			if(awgm_params.verbose){
				std::cerr << "          Bucket " << i << ": L2-norm " << r_bucket.bucket_ell2norms[i] << std::endl;
			}

			if (P_Lambda_ResNE_square >= awgm_params.alpha*resNE_norm*awgm_params.alpha*resNE_norm) {
				break;
			}
		}

		// Replace indizes from last iteration
		for(auto& lambda : LambdaTrial){
			resNE.insert(std::make_pair<decltype(lambda),T>(lambda,0.));
		}

		// New LambdaTrial: Add new indizes to u, complete to multitree
        IndexSet<Index2D> LambdaTrialDiff;
        for(auto& coeff : p){
			if(u.find(coeff.first) == u.end()){
				completeMultiTree(trialbasis, coeff.first, u, LambdaTrialDiff, 0, true);
			}
		}

        // Compute new u_leafs
        u_leafs.clear();
        FillWithZeros(LambdaTrialDiff,u_leafs);

        LambdaTrial = supp(u);

        // New LambdaTest: stable expansion of LambdaTrial (on dummy vector Ap)
        Ap.clear();
        getStableExpansion(trialbasis, testbasis, LambdaTrial, Ap);
        LambdaTest = supp(Ap);

        if(awgm_params.verbose){
            std::cout << "       LambdaTrial:  raw extension  " <<  std::setw(8) << std::right <<  p.size() << std::endl;
            std::cout << "                     multitree ext. " <<  std::setw(8) <<  LambdaTrialDiff.size() << std::endl;
            std::cout << "                     total size     " <<  std::setw(8) <<  LambdaTrial.size() << std::endl;
            std::cout << "       LambdaTest:   total size     " <<  std::setw(8) <<  LambdaTest.size() << std::left << std::endl;
            if(awgm_params.verbose_extra){
            	std::cout << LambdaTrialDiff;
            }
            std::cout << std::endl;
        }
    }

    std::cerr << "AWGM reached maximal iteration number " << awgm_params.max_its << ": "
    		  << "Residual NE = " << resNE_norm << " "
    		  << ", Residual Au-f = " << res_norm << " "
              << ", awgm_tol = " << awgm_params.tol << std::endl << std::endl;

    if(awgm_params.print_info){
    	if(awgm_params.verbose){
            std::cout << "=====>  Writing information to file " << std::endl<< std::endl;
    	}
    	awgm_info.print(awgm_params.info_filename.c_str());
    }
    if(awgm_params.plot_solution && flens::IsSame<Index,Index2D>::value){
    	if(awgm_params.verbose){
            std::cout << "=====>  Plotting solution to file " << std::endl << std::endl;
    	}
        plot2D<T,TrialBasis,TrialPrec>(trialbasis, u, trialPrec, exact_sol, 0., 1., 0., 1., 0.01, awgm_params.plot_filename.c_str());
    }
}

template <typename Index, typename TrialBasis, typename TestBasis,
		  typename LocalOperator, typename LocalOperatorTransp, typename RHS,
		  typename TrialPrec, typename TestPrec>
void
MultiTreeAWGM_PG<Index,TrialBasis,TestBasis,LocalOperator,LocalOperatorTransp,RHS,TrialPrec,TestPrec>::
set_sol(sol_fct_2d _sol)
{
	exact_sol = _sol;
}

} // namespace lawa
