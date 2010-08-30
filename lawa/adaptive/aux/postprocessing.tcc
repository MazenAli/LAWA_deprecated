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
 #include <iomanip>
 #include <lawa/righthandsides/righthandsides.h>
 #include <lawa/operators/operators.h>
 
namespace lawa {

template <typename T, typename Index, typename MA, typename RHS>
T
estimateError_Au_M_f(MA &A, RHS &F, const Coefficients<Lexicographical,T,Index> & u,
                     const IndexSet<Index> &LambdaCol)
{
    Coefficients<Lexicographical,T,Index> Au(u.d,u.d_), f(u.d,u.d_), res(u.d,u.d_);
    std::cout << "PostProcessing: Size of LambdaCol: " << LambdaCol.size() << std::endl;
    Au = mv(LambdaCol,A,u);
    std::cout << "PostProcessing: A*u finished." << std::endl;
    f  = F(LambdaCol);
    std::cout << "F finished." << std::endl;
    res = Au-f;
    std::cout << "res finished." << std::endl;
    return res.norm(2.)/f.norm(2.);
}

template <typename T, typename Index, typename MA, typename RHS>
T
estimateError_H_energy(MA &A_H, RHS &F_H, const Coefficients<Lexicographical,T,Index> & u,
                       T HNormOfExactSolution)
{
    Coefficients<Lexicographical,T,Index> Au(u.d,u.d_), f_H(u.d,u.d_);
    Au = mv(supp(u),A_H,u);
    T uAu = u*Au;
    f_H   = F_H(supp(u));
    T fu  = f_H*u;

    return std::sqrt(fabs(std::pow(HNormOfExactSolution,2)- 2*fu + uAu));

}

template<typename T, typename Preconditioner>
T
estimate_SpaceTimeError_L0T_H1(Coefficients<Lexicographical,T,Index2D> & u, 
                               Coefficients<Lexicographical,T,Index2D> & u_exact,
                               const Preconditioner &P)
{
    T error_est = 0.0;
    IndexSet<Index2D> Lambda = supp(u);
    IndexSet<Index2D> ExpandedLambda = supp(u_exact);
    
    typedef typename IndexSet<Index2D>::const_iterator const_set_it;
	
	for (const_set_it it=ExpandedLambda.begin(); it!=ExpandedLambda.end(); ++it) {
	    T prec = P((*it));
        if( Lambda.count(*it) > 0){
            //std::cout << (*it) << " u - u_exact = " << (u[*it] - u_exact[*it]) << std::endl;
            error_est += pow2i<T>(2*(*it).index2.j) * prec * (u[*it] - u_exact[*it]) * prec * (u[*it] - u_exact[*it]);
        }
        else{
            //std::cout << (*it) << " u_exact = " << u_exact[*it] << std::endl;
            error_est += pow2i<T>(2*(*it).index2.j) * prec * u_exact[*it] * prec * u_exact[*it];
        }
    }
    
    return std::sqrt(error_est);
}

template<typename T, typename Preconditioner>
T
estimate_SpaceTimeError_W0T(Coefficients<Lexicographical,T,Index2D> & u, 
                            Coefficients<Lexicographical,T,Index2D> & u_exact,
                            const Preconditioner &P)
{
    T error_est = 0.0;
    IndexSet<Index2D> Lambda = supp(u);
    IndexSet<Index2D> ExpandedLambda = supp(u_exact);
    
    typedef typename IndexSet<Index2D>::const_iterator const_set_it;
	
	for (const_set_it it=ExpandedLambda.begin(); it!=ExpandedLambda.end(); ++it) {
	    T prec = P((*it));
        if( Lambda.count(*it) > 0){
            //std::cout << (*it) << " u - u_exact = " << (u[*it] - u_exact[*it]) << std::endl;
            error_est += (pow2i<T>(2*(*it).index1.j - 2*(*it).index2.j) + pow2i<T>(2*(*it).index2.j)) 
                            * prec * (u[*it] - u_exact[*it]) * prec * (u[*it] - u_exact[*it]);
        }
        else{
            //std::cout << (*it) << " u_exact = " << u_exact[*it] << std::endl;
            error_est += (pow2i<T>(2*(*it).index1.j - 2*(*it).index2.j) + pow2i<T>(2*(*it).index2.j)) 
                            * prec * u_exact[*it] * prec * u_exact[*it];
        }
    }
    
    return std::sqrt(error_est);
}


/*template <typename T, typename Index, typename MA, typename RHS>
T
estimateError_H_energy(T t, MA &A_H, RHS &F_H, const Coefficients<Lexicographical,T,Index> & u,
					   T HNormOfExactSolution)
{
	Coefficients<Lexicographical,T,Index> Au(u.d,u.d_), f_H(u.d,u.d_);
	Au = mv(t,supp(u),A_H,u);
	T uAu = u*Au;
	f_H   = F_H(t, supp(u));
	T fu  = f_H*u;

	return std::sqrt(fabs(std::pow(HNormOfExactSolution,2)- 2*fu + uAu));

}

template <typename T, typename Index, typename MA, typename RHS>
T
estimateError_H_energy_Mario(MA &A_H, RHS &F_H, const Coefficients<Lexicographical,T,Index> & u,
					        T HNormOfExactSolution)
{
	Coefficients<Lexicographical,T,Index> Au(u.d,u.d_), f_H(u.d, u.d_), Af_H(u.d,u.d_);	
	
	Au = mv(supp(u),A_H,u);
	T uAAu = Au * Au;
	f_H   = F_H(supp(u));
	T Afu  = Au * f_H;
	
	std::cerr.precision(20);
    std::cerr << "(0.5 * uAAu - Afu) = " <<std::setw(15) <<  (0.5 * uAAu - Afu) << std::endl;

	return 2*std::sqrt(fabs(HNormOfExactSolution - (0.5 * uAAu - Afu)));

}*/

/*template<typename T, typename BasisX>
T
estimate_SpaceTimeError_L0T_H1_energy(Coefficients<Lexicographical,T,Index2D> & u,
                            BasisX& basis, T (*exact_t)(T), flens::DenseVector<flens::Array<T> > singpts_x,
                            T (*dd_exact_x)(T), flens::DenseVector<flens::Array<T> > singpts_t,
                            T (*H1_t_NormOfExactSolution)(T), T deltaT)
{
    // BilinearForm
    typedef HelmholtzOperator1D<T, BasisX> H1NormOp;
    typedef DiagonalMatrixPreconditioner1D<T, BasisX, H1NormOp> DiagPrec;  
    typedef Compression<T,BasisX, H1NormOp> PDECompression2D;	//not yet implemented
    typedef MapMatrix<T,Index1D,H1NormOp,PDECompression2D,DiagPrec> MA;

    H1NormOp a(basis, 1.);
    DiagPrec P(a);
    MA A(a, P);
    
    // right hand side
    typedef TimedepSeparableRHS1D<T,BasisX>     RhsIntegral;
    typedef RHS<T,Index1D, RhsIntegral, DiagPrec> Rhs;
    
    SeparableFunction2D<T> SepFunc(exact_t, singpts_t,
							       dd_exact_x, singpts_x);

	RhsIntegral RHS(basis, SepFunc, 10);
	Rhs F(RHS,P);

	typedef typename Coefficients<Lexicographical,T,Index2D>::const_iterator const_it;
    typedef typename Coefficients<Lexicographical,T,Index1D>::value_type val_type;
	
    T error_est = 0.0;
    // t = 0
    Coefficients<Lexicographical, T, Index1D> u_t;
    for(const_it it1 = u.begin(); it1 != u.end(); ++it1){
        T val = 0.0;
        for(const_it it2 = u.begin(); it2 != u.end(); ++it2){
            if( (*it2).first.index2 == (*it1).first.index2 ){
                if( (*it2).first.index1.xtype == XBSpline ){
                    val += (*it2).second * basis.mra.phi(0., (*it2).first.index1.j, (*it2).first.index1.k);                    
                }
                else{
                    val += (*it2).second * basis.psi(0., (*it2).first.index1.j, (*it2).first.index1.k);                                        
                }
            }
        }
        u_t.insert(val_type((*it1).first.index2, val));
    }
    error_est += estimateError_H_energy(0., A, F, u_t, H1_t_NormOfExactSolution(0.));
    
    for(T t = deltaT; t < 1.; t += deltaT){
        u_t.erase(u_t.begin(), u_t.end());
        for(const_it it1 = u.begin(); it1 != u.end(); ++it1){
            T val = 0.0;
            for(const_it it2 = u.begin(); it2 != u.end(); ++it2){
                if( (*it2).first.index2 == (*it1).first.index2 ){
                    if( (*it2).first.index1.xtype == XBSpline ){
                        val += (*it2).second * basis.mra.phi(t, (*it2).first.index1.j, (*it2).first.index1.k);                    
                    }
                    else{
                        val += (*it2).second * basis.psi(t, (*it2).first.index1.j, (*it2).first.index1.k);                                        
                    }
                }
            }
            u_t.insert(val_type((*it1).first.index2, val));
        }
        error_est += estimateError_H_energy(t, A, F, u_t, H1_t_NormOfExactSolution(t));
    }
    
    // t = 1
    u_t.erase(u_t.begin(), u_t.end());
    for(const_it it1 = u.begin(); it1 != u.end(); ++it1){
        T val = 0.0;
        for(const_it it2 = u.begin(); it2 != u.end(); ++it2){
            if( (*it2).first.index2 == (*it1).first.index2 ){
                if( (*it2).first.index1.xtype == XBSpline ){
                    val += (*it2).second * basis.mra.phi(1., (*it2).first.index1.j, (*it2).first.index1.k);                    
                }
                else{
                    val += (*it2).second * basis.psi(1., (*it2).first.index1.j, (*it2).first.index1.k);                                        
                }
            }
        }
        u_t.insert(val_type((*it1).first.index2, val));
    }
    error_est += estimateError_H_energy(1., A, F, u_t, H1_t_NormOfExactSolution(1.));
    
    return std::sqrt(error_est);
}*/

} // namespace lawa
