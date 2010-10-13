/*
  LAWA - Library for Adaptive Wavelet Applications.
  Copyright (C) 2008,2009  Sebastian Kestler, Mario Rometsch, Kristina Steih, Alexander Stippler.

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
#include <lawa/adaptive/aux/timer.h>

namespace lawa {


template <typename T, typename Index, typename Basis, typename MA, typename RHS>
S_ADWAV<T,Index,Basis,MA,RHS>::S_ADWAV(const Basis &_basis, MA &_A, RHS &_F, T _contraction,
                                 T start_threshTol, T start_linTol, T start_resTol,
                                 int _NumOfIterations, int _MaxItsPerThreshTol, T _eps)
    : basis(_basis), A(_A), F(_F), contraction(_contraction), threshTol(start_threshTol), linTol(start_linTol),
      resTol(start_resTol), NumOfIterations(_NumOfIterations), MaxItsPerThreshTol(_MaxItsPerThreshTol), eps(_eps)
{
    solutions.resize(NumOfIterations);
    residuals.resize(NumOfIterations);
    times.resize(NumOfIterations);
}

template <typename T, typename Index, typename Basis, typename MA, typename RHS>
void
S_ADWAV<T,Index,Basis,MA,RHS>::solve_cg(const IndexSet<Index> &InitialLambda, T H1norm)
{
    Timer timer;
    
    int d=InitialLambda.d, d_=InitialLambda.d_;
    IndexSet<Index> LambdaActive(d,d_), LambdaThresh(d,d_), LambdaActivable(d,d_), DeltaLambda(d,d_);
    Coefficients<Lexicographical,T, Index> u(d,d_), f(d,d_), Au(d,d_), r(d,d_);

    LambdaActive = InitialLambda;
    T old_res = 0.;
    int its_per_threshTol=0;
    std::cout << "Simple adaptive solver started." << std::endl;
    std::stringstream filename;
    filename << "s-adwav-realline-helmholtz-otf_" << d << "_" << d_ << ".dat";
    std::ofstream file(filename.str().c_str());

    for (int its=0; its<NumOfIterations; ++its) {
    	std::cout << "*** " << its+1 << ".iteration" << std::endl;
        
        timer.start();
        
        //Initialization step
        FillWithZeros(LambdaActive,u);
        f = F(LambdaActive);
        T f_norm_LambdaActive = f.norm(2.);

        //Galerkin step
        T r_norm_LambdaActive = 0.0;
        std::cout << "   CG solver started with N = " << LambdaActive.size() << std::endl;
        int iterations = CG_Solve(LambdaActive, A, u, f, r_norm_LambdaActive, linTol);
        std::cout << "   ...finished." << std::endl;

        //Threshold step
        u = THRESH(u,threshTol);
        solutions[its] = u;
        LambdaThresh = supp(u);
        std::cout << "    Size of thresholded u = " << LambdaThresh.size() << std::endl;

        timer.stop();
        T time1 = timer.elapsed();
        T Error_H_energy = 0.;
        if (H1norm>0) Error_H_energy = estimateError_H_energy(A, F, u, H1norm);

        timer.start();
        //Computing residual
        DeltaLambda = C(LambdaThresh, contraction, basis);
        std::cout << "   Computing rhs for DeltaLambda (size = " << DeltaLambda.size() << ")" << std::endl;
        f = F(DeltaLambda);
        std::cout << "   ...finished" << std::endl;
        T f_norm_DeltaLambda = f.norm(2.);
        std::cout << "   Computing residual for DeltaLambda (size = " << DeltaLambda.size() << ")" << std::endl;
        //Au = mv(DeltaLambda,A,u);
        Au = mv_sparse(DeltaLambda,A,u);
        r  = Au-f;
        T r_norm_DeltaLambda = r.norm(2.);
        T numerator   = r_norm_DeltaLambda*r_norm_DeltaLambda + r_norm_LambdaActive*r_norm_LambdaActive;
        T denominator = f_norm_DeltaLambda*f_norm_DeltaLambda + f_norm_LambdaActive*f_norm_LambdaActive;
        T estim_res   = std::sqrt(numerator/denominator);
        std::cout << "   ...finished" << std::endl;
        residuals[its] = estim_res;

        file << LambdaThresh.size() << " " << estim_res << " " << Error_H_energy << std::endl;

        r = THRESH(r,threshTol);
        LambdaActive = LambdaThresh+supp(r);

        //Check if residual is decreasing, if not decrease threshold tolerance
        if (fabs(estim_res-old_res)<resTol || its_per_threshTol>MaxItsPerThreshTol) {
            threshTol *= 0.5;
            linTol      *= 0.5;
            resTol    *= 0.5;
            its_per_threshTol = 0;
        }
        ++its_per_threshTol;
        old_res = estim_res;
        timer.stop();
        if (its==0) times[its] = time1+timer.elapsed();
        else		times[its] = times[its-1] + time1 + timer.elapsed();
        
        std::cout << "S-ADWAV: " << its+1 << ".iteration: Size of Lambda = " << supp(u).size() << ", cg-its = " << iterations
                  << ", residual = " << estim_res << " , current threshTol = " << threshTol << std::endl << std::endl;

    }
}

template <typename T, typename Index, typename Basis, typename MA, typename RHS>
void
S_ADWAV<T,Index,Basis,MA,RHS>::solve_gmres(const IndexSet<Index> &InitialLambda)
{
    Timer timer;
    
	int d=InitialLambda.d, d_=InitialLambda.d_;
	IndexSet<Index> LambdaActive(d,d_), LambdaThresh(d,d_), LambdaActivable(d,d_), DeltaLambda(d,d_);
	Coefficients<Lexicographical,T, Index> u(d,d_), f(d,d_), Au(d,d_), r(d,d_);

	LambdaActive = InitialLambda;
	T old_res = 0.;
	int its_per_threshTol=0;
	std::cout << "Simple adaptive solver started." << std::endl;
	for (int its=0; its<NumOfIterations; ++its) {
       
        timer.start();
	
		//Initialization step
		FillWithZeros(LambdaActive,u);
		f = F(LambdaActive);
		T f_norm_LambdaActive = f.norm(2.);

		//Galerkin step
		T r_norm_LambdaActive = 0.0;
		std::cout << "   GMRES solver started with N = " << LambdaActive.size() << std::endl;
		int iterations = GMRES_Solve(LambdaActive, A, u, f, r_norm_LambdaActive, linTol);
		std::cout << "   ...finished." << std::endl;

		//Threshold step
		u = THRESH(u,threshTol);
		solutions[its] = u;
		LambdaThresh = supp(u);

		//Computing residual
		DeltaLambda = C(LambdaThresh, contraction, basis);
		std::cout << "   Computing rhs for DeltaLambda (size = " << DeltaLambda.size() << ")" << std::endl;
		f = F(DeltaLambda);
		std::cout << "   ...finished" << std::endl;
		T f_norm_DeltaLambda = f.norm(2.);
		std::cout << "   Computing residual for DeltaLambda (size = " << DeltaLambda.size() << ")" << std::endl;
		Au = mv_sparse(DeltaLambda,A,u);
		r  = Au-f;
		T r_norm_DeltaLambda = r.norm(2.);
		T numerator   = r_norm_DeltaLambda*r_norm_DeltaLambda + r_norm_LambdaActive*r_norm_LambdaActive;
		T denominator = f_norm_DeltaLambda*f_norm_DeltaLambda + f_norm_LambdaActive*f_norm_LambdaActive;
		T estim_res   = std::sqrt(numerator/denominator);
		std::cout << "   ...finished" << std::endl;
		residuals[its] = estim_res;

		r = THRESH(r,threshTol);
        std::cout << "r_Thres: " << r << std::endl;
		//LambdaActive = LambdaThresh+supp(r);
		LambdaActive = LambdaActive+supp(r);

		//Check if residual is decreasing, if not decrease threshold tolerance
		if (fabs(estim_res-old_res)<resTol || its_per_threshTol>MaxItsPerThreshTol) {
		//if(old_res - estim_res < resTol){
			threshTol *= 0.5;
			linTol	  *= 0.5;
			resTol    *= 0.5;
			its_per_threshTol = 0;
		}
		++its_per_threshTol;
		old_res = estim_res;
        timer.stop();
        times[its] = timer.elapsed();
        
		std::cout << "S-ADWAV: " << its+1 << ".iteration: Size of Lambda = " << supp(u).size() << ", gmres-its = " << iterations;
		std::cout << ", residual = " << estim_res << " , current threshTol = " << threshTol << std::endl;
	}
}

template <typename T, typename Index, typename Basis, typename MA, typename RHS>
void
S_ADWAV<T,Index,Basis,MA,RHS>::solve_cgls(const IndexSet<Index> &InitialLambda)
{
    Timer timer;

    int d=InitialLambda.d, d_=InitialLambda.d_;
    IndexSet<Index> LambdaActive(d,d_), LambdaThresh(d,d_), LambdaActivable(d,d_);
    IndexSet<Index> 	 LambdaActive_test(d,d_), DeltaLambda(d,d_);
    Coefficients<Lexicographical,T, Index>      u(d,d_), f(d,d_), Au(d,d_), r(d,d_);

    LambdaActive = InitialLambda;
    T old_res = 0.;
    int its_per_threshTol=0;
    std::cout << "Simple adaptive cgls time solver started." << std::endl;
    std::stringstream filename;
    filename << "s-adwav-cgls_" << d << "_" << d_ << ".dat";
    std::ofstream file(filename.str().c_str());

    for (int its=0; its<NumOfIterations; ++its) {
    	std::cout << "*** " << its+1 << ".iteration" << std::endl;

        timer.start();

        //Initialization step
        LambdaActive_test = LambdaActive;// + C(LambdaActive,contraction,basis);
        FillWithZeros(LambdaActive,u);
        f  = F (LambdaActive_test);

        //Galerkin step
        T r_norm_LambdaActive = 0.0;
        std::cout << "   CG solver started with N = " << LambdaActive.size() << std::endl;
        int iterations = CGLS_Solve(LambdaActive_test, LambdaActive,
									A, u, f, r_norm_LambdaActive, linTol, 100000);
        std::cout << "   ...finished after " << iterations << " iterations with residual = " << r_norm_LambdaActive << std::endl;

        //Threshold step
        u = THRESH(u,threshTol);
        solutions[its] = u;
        LambdaThresh = supp(u);
        std::cout << "    Size of thresholded u = " << LambdaThresh.size() << std::endl;

        //Computing residual for operator part
        DeltaLambda = LambdaThresh + C(LambdaThresh, contraction, basis);
        std::cout << "   Computing rhs for DeltaLambda_test (size = " << DeltaLambda.size() << ")" << std::endl;
        f = F(DeltaLambda);
        std::cout << "   ...finished" << std::endl;
        T f_norm_DeltaLambda = f.norm(2.);
        std::cout << "   Computing residual for DeltaLambda_test_operator..." << std::endl;
        Au = mv_sparse(DeltaLambda,A,u);
        r  = Au-f;
        T r_norm_DeltaLambda = r.norm(2.);
        std::cout << "   ...finished" << std::endl;

        T numerator   = r_norm_DeltaLambda*r_norm_DeltaLambda;
        T denominator = f_norm_DeltaLambda*f_norm_DeltaLambda;
        T estim_res   = std::sqrt(numerator/denominator);
        std::cout << "   ...finished" << std::endl;
        residuals[its] = estim_res;

        file << u.size() << " " << LambdaActive.size() << " " << estim_res << " " << iterations << " " << r_norm_LambdaActive << std::endl;

        r = THRESH(r,threshTol);
        LambdaActive = LambdaThresh+supp(r);

        //Check if residual is decreasing, if not decrease threshold tolerance
        if (fabs(estim_res-old_res)<resTol || its_per_threshTol>MaxItsPerThreshTol) {
            threshTol *= 0.5;
            linTol      *= 0.5;
            //resTol    *= 0.5;
            its_per_threshTol = 0;
        }
        ++its_per_threshTol;
        old_res = estim_res;
        timer.stop();
        if (its==0) times[its] = timer.elapsed();
        else		times[its] = times[its-1] + timer.elapsed();

        std::cout << "S-ADWAV: " << its+1 << ".iteration: Size of Lambda = " << supp(u).size() << ", cgls-its = " << iterations
                  << ", residual = " << estim_res << " , current threshTol = " << threshTol << std::endl << std::endl;

    }
}

}    //namespace lawa

