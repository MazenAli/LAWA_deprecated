#include <lawa/methods/adaptive/algorithms/algorithms.h>

#include <ios>
#include <iostream>
#include <fstream>
#include <iomanip>

namespace lawa {

AWGM_Parameters::AWGM_Parameters(double _tol, double _alpha, size_t _max_its,
		size_t _max_basissize, bool _reset_resNE, bool _reset_res, bool _print_info,
		bool _verbose, bool _plot_solution, bool _verbose_extra, size_t _hashmapsize)
: tol(_tol), alpha(_alpha), max_its(_max_its), max_basissize(_max_basissize),
  reset_resNE(_reset_resNE), reset_res(_reset_res), print_info(_print_info),
  verbose(_verbose), plot_solution(_plot_solution), verbose_extra(_verbose_extra),
  hashmapsize(_hashmapsize)
{}

void
AWGM_Parameters::print()
{
	std::cout << "###### AWGM Parameters #################" << std::endl;
	std::cout << std::left << std::setw(18) << "# tol:" << std::setw(16) <<  tol << std::endl;
	std::cout << std::left << std::setw(18) << "# alpha:" << std::setw(16) <<  alpha << std::endl;
	std::cout << std::left << std::setw(18) << "# max_its:" << std::setw(16) <<  max_its << std::endl;
	std::cout << std::left << std::setw(18) << "# max_basissize:" << std::setw(16) <<  max_basissize << std::endl;
	std::cout << std::left << std::setw(18) << "# reset_resNE:" << std::setw(16) <<  (reset_resNE?"true":"false") << std::endl;
	std::cout << std::left << std::setw(18) << "# reset_res:" << std::setw(16) <<  (reset_res?"true":"false") << std::endl;
	std::cout << std::left << std::setw(18) << "# print_info:" << std::setw(16) <<  (print_info?"true":"false") << std::endl;
	std::cout << std::left << std::setw(18) << "# verbose:" << std::setw(16) <<  (verbose?"true":"false") << (verbose_extra?" (extra)":"") << std::endl;
	std::cout << std::left << std::setw(18) << "# plot_solution:" << std::setw(16) <<  (plot_solution?"true":"false") << std::endl;
	std::cout << std::left << std::setw(18) << "# hashmapsize:" << std::setw(16) <<  hashmapsize << std::endl;
	std::cout << "#########################################" << std::endl << std::endl;
}


IS_Parameters::IS_Parameters(bool _adaptive_tol, size_t _max_its, double _init_tol,
			  double _res_reduction, double _absolute_tol, bool _verbose)
: adaptive_tol(_adaptive_tol), max_its(_max_its), init_tol(_init_tol),
  res_reduction(_res_reduction), absolute_tol(_absolute_tol), verbose(_verbose)
{}

void
IS_Parameters::print()
{
	std::cout << "###### Inner Solver Parameters ##########" << std::endl;
	std::cout << std::left << std::setw(18) << "# adaptive_tol:" << std::setw(16) <<  (adaptive_tol?"true":"false") << std::endl;
	std::cout << std::left << std::setw(18) << "# max_its:" << std::setw(16) <<  max_its << std::endl;
	std::cout << std::left << std::setw(18) << "# init_tol:" << std::setw(16) <<  init_tol << std::endl;
	std::cout << std::left << std::setw(18) << "# res_reduction:" << std::setw(16) <<  res_reduction << std::endl;
	std::cout << std::left << std::setw(18) << "# absolute_tol:" << std::setw(16) <<  absolute_tol << std::endl;
	std::cout << std::left << std::setw(18) << "# verbose:" << std::setw(16) <<  (verbose?"true":"false") << std::endl;
	std::cout << "#########################################" << std::endl << std::endl;
}

void
AWGM_Information::print(const char* filename)
{
    std::ofstream infofile(filename);
    if(infofile.is_open()){
    	infofile << "# It Res Res_NE SizeTrial SizeTest SizeTestResNE SizeTestRes" << std::endl;
    	for(size_t i=0; i < awgm_res.size(); ++i){
    		infofile << i << " " << awgm_res[i] << " " << awgm_resNE[i] << " "
    		        << sizeLambdaTrial[i] << " " << sizeLambdaTest[i] << " "
    		        << sizeLambdaResNE[i] << " " << sizeLambdaRes[i] << std::endl;
    	}
        infofile.close();
    }
    else{
    	std::cerr << "Error opening file " << filename << " for writing! " << std::endl;
    }
}

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
				AWGM_Parameters& _awgm_params, IS_Parameters& _is_params)
: awgm_params(_awgm_params), is_params(_is_params),
  trialbasis(_trialbasis), testbasis(_testbasis), Op(_Op), OpTransp(_OpTransp), F(_F),
  trialPrec(_trialPrec), testPrec(_testPrec), exact_sol(nullptr)
{}

template <typename Index, typename TrialBasis, typename TestBasis,
		  typename LocalOperator, typename LocalOperatorTransp, typename RHS,
		  typename TrialPrec, typename TestPrec>
void
MultiTreeAWGM_PG<Index,TrialBasis,TestBasis,LocalOperator,LocalOperatorTransp,RHS,TrialPrec,TestPrec>::
cgls_solve(Coefficients<Lexicographical,T,Index> &u)
{
    IndexSet<Index> LambdaTrial, LambdaTest;
	if(u.size() > 0){
		LambdaTrial = supp(u);
	}
	else{
		if(flens::IsSame<Index,Index2D>::value){
			getSparseGridIndexSet(trialbasis,LambdaTrial,2,0);
		}
		else{
			std::cerr << "Need some initial set in trial space! " << std::endl;
			exit(1);
		}
	}
    Coefficients<Lexicographical,T,Index2D> Lambda_aux;
	getStableExpansion(trialbasis, testbasis, LambdaTrial, Lambda_aux);
	LambdaTest = supp(Lambda_aux);

	cgls_solve(u, LambdaTrial, LambdaTest);
}

template <typename Index, typename TrialBasis, typename TestBasis,
		  typename LocalOperator, typename LocalOperatorTransp, typename RHS,
		  typename TrialPrec, typename TestPrec>
void
MultiTreeAWGM_PG<Index,TrialBasis,TestBasis,LocalOperator,LocalOperatorTransp,RHS,TrialPrec,TestPrec>::
cgls_solve(Coefficients<Lexicographical,T,Index> &u, IndexSet<Index>& LambdaTrial, IndexSet<Index>& LambdaTest)
{
    //---------------------------------------//
    //------- AWGM Initialization -----------//
    //---------------------------------------//

    /// Initialization of hash map vectors
    Coefficients<Lexicographical,T,Index2D> f(SIZEHASHINDEX2D),
    										r(SIZEHASHINDEX2D),
    										s(SIZEHASHINDEX2D),
                                            p(SIZEHASHINDEX2D),
                                            q(SIZEHASHINDEX2D),
                                            Ap(SIZEHASHINDEX2D),
                                            res(SIZEHASHINDEX2D),     // approximate residual for f-Au
                                    		resNE(SIZEHASHINDEX2D),   // cone around u_leafs
                                            u_leafs(SIZEHASHINDEX2D); // "leafs" of u

	FillWithZeros(LambdaTrial,u_leafs);
	FillWithZeros(LambdaTest,res);

    // Use coefficient vectors for preconditioners
    // (no possibility to store precomputed values in the actual precs yet!)
    Coefficients<Lexicographical,T,Index2D> leftP(SIZEHASHINDEX2D);
    Coefficients<Lexicographical,T,Index2D> rightP(SIZEHASHINDEX2D);

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

        // Calculate right hand side
		//f = F(LambdaTest);

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

		// Compute Preconditioner values
		for(auto& lambda : LambdaTest){
			if(leftP.find(lambda)==leftP.end()){
				leftP.insert(std::make_pair<decltype(lambda),T>(lambda, testPrec(lambda)));
			}
		}
		for(auto& lambda : LambdaTrial){
			if(rightP.find(lambda)==rightP.end()){
				rightP.insert(std::make_pair<decltype(lambda),T>(lambda, trialPrec(lambda)));
			}
		}

		// Initial step
	    Op.eval(u,r,rightP,leftP);
		//r -= f;
	    for(auto& lambda : LambdaTest){
	    	r[lambda] -= leftP[lambda]*F(lambda);
	    }
		r *= -1;
		OpTransp.eval(r,s,leftP,rightP);
		p = s;
		gamma_cgls_Prev = s*s;

		// CGLS Iterations
		for(size_t cgls_its=0; cgls_its <= is_params.max_its; ++cgls_its){

			Ap.setToZero();						// Ap = A*p
			Op.eval(p,Ap,rightP,leftP);

			alpha = gamma_cgls_Prev / (Ap*Ap);
			u += alpha * p;
			r -= alpha * Ap;

			s.setToZero();
			OpTransp.eval(r,s,leftP,rightP);	// s = A^T*r

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
				break;
			}

			beta  = gamma_cgls/gamma_cgls_Prev;
			p *= beta;
			p += s;
			gamma_cgls_Prev = gamma_cgls;
		}

        //---------------------------------------//
        //---- COMPUTE APPROXIMATE RESIDUAL -----//
        //---------------------------------------//

		// Compute cone around solution u
		IndexSet<Index2D>	Cdiff_u_leafs;	// new indizes in cone
        res.setToZero();
        if(awgm_params.reset_resNE){
        	resNE.clear();
        }
        extendMultiTree(trialbasis, u, resNE, Cdiff_u_leafs, "standard", false, true);

        // Compute stable expansion in test basis
        IndexSet<Index2D> LambdaResNE = supp(resNE);
        if(awgm_params.reset_res){
        	res.clear();
        }
        getStableExpansion(trialbasis, testbasis, LambdaResNE, res);

        // Compute rhs on expanded test index set
        IndexSet<Index2D> LambdaRes = supp(res);
        //f = F(LambdaRes);

        // Update Preconditioner values
		for(auto& lambda : LambdaRes){
			if(leftP.find(lambda)==leftP.end()){
				leftP.insert(std::make_pair<decltype(lambda),T>(lambda, testPrec(lambda)));
			}
		}
		for(auto& lambda : LambdaResNE){
			if(rightP.find(lambda)==rightP.end()){
				rightP.insert(std::make_pair<decltype(lambda),T>(lambda, trialPrec(lambda)));
			}
		}

		// Compute residual Au - f on expanded test index set
        res.setToZero();
        Op.eval(u,res,rightP,leftP); 	// res = A*u
        //res -= f;
	    for(auto& lambda : LambdaRes){
	    	res[lambda] -= leftP[lambda]*F(lambda);
	    }
        res_norm = res.norm(2.);

        awgm_info.awgm_res.push_back(res_norm);
        awgm_info.sizeLambdaRes.push_back(LambdaRes.size());

        // Compute residual of NE: A^TA*u - A^T*f
        resNE.setToZero();
        OpTransp.eval(res,resNE,leftP,rightP); 	// resNE = A^T*res
        resNE_norm = resNE.norm(2.);

        awgm_info.awgm_resNE.push_back(resNE_norm);
        awgm_info.sizeLambdaResNE.push_back(LambdaResNE.size());

        if(awgm_params.verbose){
            std::cout << "   --- Approximate Residual ---" << std::endl;
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
                      << ", awgm_tol = " << awgm_params.tol << std::endl;

            if(awgm_params.print_info){
            	awgm_info.print();
            }

            return;
        }

        //---------------------------------------//
        //---- COMPUTE NEXT INDEX SET -----------//
        //---------------------------------------//

        // Remove indizes from last iteration
		for(auto& lambda : LambdaTrial){
			resNE.erase(lambda);
		}

		// Compute buckets
		T threshbound = std::sqrt(1.-awgm_params.alpha*awgm_params.alpha) * resNE_norm/std::sqrt(T(resNE.size()));
		Coefficients<Bucket,T,Index2D> r_bucket;
		r_bucket.bucketsort(resNE, threshbound);

		// Add buckets to dummy vector: p
		T P_Lambda_ResNE_square = gamma_cgls;

        if(awgm_params.verbose){
            std::cout << "   --- Computing next index set ---" << std::endl;
            std::cout << "       Threshbound:         " << threshbound <<  std::endl;
            std::cout << "       ||P_{Lambda}r ||_2:  " << std::sqrt(P_Lambda_ResNE_square) << std::endl;
            std::cout << "       alpha*Residual_NE:   " << awgm_params.alpha*resNE_norm << std::endl;
        }

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
            	std::cout << LambdaTrialDiff << std::endl << std::endl;
            }
        }
    }

    std::cerr << "AWGM reached maximal iteration number " << awgm_params.max_its << ": "
    		  << "Residual NE = " << resNE_norm << " "
    		  << ", Residual Au-f = " << res_norm << " "
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
        plot2D<T,TrialBasis,TrialPrec>(trialbasis, u, trialPrec, exact_sol, 0., 1., 0., 1., 0.01, "awgm_cgls_u_plot");
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
