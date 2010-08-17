namespace lawa{

    
// OPERATOR_LHSMATRIX
template<typename T, typename Basis, typename Problem, typename BilinearForm, typename RHSIntegral>
ThetaScheme<T, Basis, Problem, BilinearForm, RHSIntegral>::Operator_LHSMatrix::
Operator_LHSMatrix(ThetaScheme<T, Basis, Problem, BilinearForm, RHSIntegral>* _scheme, 
                  const Basis& _basis, const BilinearForm& _a)
    : basis(_basis), a(_a)
{   
    scheme = _scheme;
}

template<typename T, typename Basis, typename Problem, typename BilinearForm, typename RHSIntegral>
T 
ThetaScheme<T, Basis, Problem, BilinearForm, RHSIntegral>::Operator_LHSMatrix::
operator()(XType xtype1, int j1, int k1,
           XType xtype2, int j2, int k2) const
{
    // (M + deltaT * theta * A_k+1)
    T timestep = time_new - time_old;
    if(xtype1 == XBSpline){
        if(xtype2 == XBSpline){
            return scheme->integral_sfsf(j1,k1,j2,k2) 
                    + timestep*scheme->theta*a(time_new, xtype1,j1,k1, xtype2, j2,k2);
        }
        else{
            return scheme->integral_sfw(j1,k1,j2,k2) 
                    + timestep*scheme->theta*a(time_new, xtype1,j1,k1, xtype2, j2,k2);            
        }
    }
    else{
        if(xtype2 == XBSpline){
            return scheme->integral_wsf(j1,k1,j2,k2) 
                    + timestep*scheme->theta*a(time_new, xtype1,j1,k1, xtype2, j2,k2);
            
        }
        else{
            return scheme->integral_ww(j1,k1,j2,k2) 
                    + timestep*scheme->theta*a(time_new, xtype1,j1,k1, xtype2, j2,k2);
        }
    }
    
}


// OPERATOR_RHSMATRIX
template<typename T, typename Basis, typename Problem, typename BilinearForm, typename RHSIntegral>
ThetaScheme<T, Basis, Problem, BilinearForm, RHSIntegral>::Operator_RHSMatrix::
Operator_RHSMatrix(const ThetaScheme<T, Basis, Problem, BilinearForm, RHSIntegral>* _scheme, 
                   const Basis& _basis, const BilinearForm& _a)
    : basis(_basis), a(_a)
{
     scheme = _scheme;
}


template<typename T, typename Basis, typename Problem, typename BilinearForm, typename RHSIntegral>
T 
ThetaScheme<T, Basis, Problem, BilinearForm, RHSIntegral>::Operator_RHSMatrix::
operator()(XType xtype1, int j1, int k1,
           XType xtype2, int j2, int k2) const
{
   // (M - deltaT * (1-theta) * A_k)
   T timestep = time_new - time_old;
   if(xtype1 == XBSpline){
       if(xtype2 == XBSpline){
           return scheme->integral_sfsf(j1,k1,j2,k2) 
                    - timestep*(1.-scheme->theta)*a(time_old, xtype1,j1,k1, xtype2, j2,k2);
       }
       else{
           return scheme->integral_sfw(j1,k1,j2,k2) 
                    - timestep*(1.-scheme->theta)*a(time_old, xtype1,j1,k1, xtype2, j2,k2);            
       }
   }
   else{
       if(xtype2 == XBSpline){
           return scheme->integral_wsf(j1,k1,j2,k2) 
                    - timestep*(1.-scheme->theta)*a(time_old, xtype1,j1,k1, xtype2, j2,k2);
       }
       else{
           return scheme->integral_ww(j1,k1,j2,k2) 
                    - timestep*(1.-scheme->theta)*a(time_old, xtype1,j1,k1, xtype2, j2,k2);
       }
   }

}

// OPERATOR_RHSVECTOR
template<typename T, typename Basis, typename Problem, typename BilinearForm, typename RHSIntegral>
ThetaScheme<T, Basis, Problem, BilinearForm, RHSIntegral>::Operator_RHSVector::
Operator_RHSVector(const ThetaScheme<T, Basis, Problem, BilinearForm, RHSIntegral>* _scheme, 
                   const Basis& _basis, const RHSIntegral& _rhs)
    : basis(_basis), rhs(_rhs)
{
     scheme = _scheme;
}

template<typename T, typename Basis, typename Problem, typename BilinearForm, typename RHSIntegral>
T 
ThetaScheme<T, Basis, Problem, BilinearForm, RHSIntegral>::Operator_RHSVector::
operator()(XType xtype, int j, int k) const
{   
    // deltaT * (theta*f_k+1 - (1-theta)*f_k)
    T timestep = time_new - time_old;
    return timestep*(scheme->theta * rhs(time_new, xtype, j, k) 
                    + (1. - scheme->theta)*rhs(time_old, xtype, j, k));
}



// THETASCHEME
template<typename T, typename Basis, typename Problem, typename BilinearForm, typename RHSIntegral>
ThetaScheme<T, Basis, Problem, BilinearForm, RHSIntegral>::
ThetaScheme(const T _theta, const Basis& _basis, const BilinearForm& _a, const RHSIntegral& _rhs)
    : theta(_theta), basis(_basis), problem(basis), phi(basis.mra), psi(basis),
      integral_sfsf(phi, phi), integral_sfw(phi, psi), integral_wsf(psi, phi), integral_ww(psi, psi),
      op_LHSMatrix(this, basis, _a), op_RHSMatrix(this, basis, _a), op_RHSVector(this, basis, _rhs)
{   
}
 
template<typename T, typename Basis, typename Problem, typename BilinearForm, typename RHSIntegral>
flens::DenseVector<flens::Array<T> > 
ThetaScheme<T, Basis, Problem, BilinearForm, RHSIntegral>::solve(T time_old, T time_new, flens::DenseVector<flens::Array<T> > u_init, int level)
{   
    op_LHSMatrix.setTimes(time_old, time_new);
    op_RHSMatrix.setTimes(time_old, time_new);
    op_RHSVector.setTimes(time_old, time_new);
    flens::SparseGeMatrix<flens::CRS<T,flens::CRS_General> > lhsmatrix = problem.getStiffnessMatrix(op_LHSMatrix, level);
    flens::SparseGeMatrix<flens::CRS<T,flens::CRS_General> > rhsmatrix = problem.getStiffnessMatrix(op_RHSMatrix, level);
    flens::DenseVector<flens::Array<T> > rhsvector = problem.getRHS(op_RHSVector, level);
    flens::DenseVector<flens::Array<T> > rhs = rhsmatrix * u_init + rhsvector;
    flens::DenseVector<flens::Array<T> > u(basis.mra.rangeI(level));
    cg(lhsmatrix, u, rhs);
    
    //std::cout << "u(" << time_new << "): " << u << std::endl; 
    return u;
} 
  
}