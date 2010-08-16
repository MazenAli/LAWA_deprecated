namespace lawa{

    
// OPERATOR_LHSMATRIX
template<typename T, typename Basis, typename Problem, typename BilinearForm, typename RHSIntegral>
ThetaScheme<T, Basis, Problem, BilinearForm, RHSIntegral>::Operator_LHSMatrix::
Operator_LHSMatrix(ThetaScheme<T, Basis, Problem, BilinearForm, RHSIntegral>* _scheme, const Basis& _basis, const BilinearForm& _a, T _timestep)
    : basis(_basis), a(_a), timestep(_timestep), time(0)
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
    if(xtype1 == XBSpline){
        if(xtype2 == XBSpline){
            return scheme->integral_sfsf(j1,k1,j2,k2) + timestep*scheme->theta*a(time, xtype1,j1,k1, xtype2, j2,k2);
        }
        else{
            return scheme->integral_sfw(j1,k1,j2,k2) + timestep*scheme->theta*a(time, xtype1,j1,k1, xtype2, j2,k2);            
        }
    }
    else{
        if(xtype2 == XBSpline){
            return scheme->integral_wsf(j1,k1,j2,k2) + timestep*scheme->theta*a(time, xtype1,j1,k1, xtype2, j2,k2);
        }
        else{
            return scheme->integral_ww(j1,k1,j2,k2) + timestep*scheme->theta*a(time, xtype1,j1,k1, xtype2, j2,k2);
        }
    }
    
}


// OPERATOR_RHSMATRIX
template<typename T, typename Basis, typename Problem, typename BilinearForm, typename RHSIntegral>
ThetaScheme<T, Basis, Problem, BilinearForm, RHSIntegral>::Operator_RHSMatrix::
Operator_RHSMatrix(const ThetaScheme<T, Basis, Problem, BilinearForm, RHSIntegral>* _scheme, const Basis& _basis, const BilinearForm& _a, T _timestep)
    : basis(_basis), a(_a), timestep(_timestep), time(0)
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
   if(xtype1 == XBSpline){
       if(xtype2 == XBSpline){
           return scheme->integral_sfsf(j1,k1,j2,k2) - timestep*(1.-scheme->theta)*a(time-timestep, xtype1,j1,k1, xtype2, j2,k2);
       }
       else{
           return scheme->integral_sfw(j1,k1,j2,k2) - timestep*(1.-scheme->theta)*a(time-timestep, xtype1,j1,k1, xtype2, j2,k2);            
       }
   }
   else{
       if(xtype2 == XBSpline){
           return scheme->integral_wsf(j1,k1,j2,k2) - timestep*(1.-scheme->theta)*a(time-timestep, xtype1,j1,k1, xtype2, j2,k2);
       }
       else{
           return scheme->integral_ww(j1,k1,j2,k2) - timestep*(1.-scheme->theta)*a(time-timestep, xtype1,j1,k1, xtype2, j2,k2);
       }
   }

}

// OPERATOR_RHSVECTOR
template<typename T, typename Basis, typename Problem, typename BilinearForm, typename RHSIntegral>
ThetaScheme<T, Basis, Problem, BilinearForm, RHSIntegral>::Operator_RHSVector::
Operator_RHSVector(const ThetaScheme<T, Basis, Problem, BilinearForm, RHSIntegral>* _scheme, const Basis& _basis, const RHSIntegral& _rhs, T _timestep)
    : basis(_basis), rhs(_rhs), timestep(_timestep), time(0) 
{
     scheme = _scheme;
}

template<typename T, typename Basis, typename Problem, typename BilinearForm, typename RHSIntegral>
T 
ThetaScheme<T, Basis, Problem, BilinearForm, RHSIntegral>::Operator_RHSVector::
operator()(XType xtype, int j, int k) const
{   
    // deltaT * (theta*f_k+1 - (1-theta)*f_k)
    return timestep*(scheme->theta * rhs(time, xtype, j, k)- (1. - scheme->theta)*rhs(time-timestep, xtype, j, k));
}



// THETASCHEME
template<typename T, typename Basis, typename Problem, typename BilinearForm, typename RHSIntegral>
ThetaScheme<T, Basis, Problem, BilinearForm, RHSIntegral>::
ThetaScheme(const T _theta, const Basis& _basis, const BilinearForm& _a, const RHSIntegral& _rhs, T _timestep)
    : basis(_basis), problem(basis), phi(basis.mra), psi(basis),
      integral_sfsf(phi, phi), integral_sfw(phi, psi), integral_wsf(psi, phi), integral_ww(psi, psi),
      op_LHSMatrix(this, basis, _a, _timestep), op_RHSMatrix(this, basis, _a, _timestep),
      op_RHSVector(this, basis, _rhs, _timestep)
{   
}
 
template<typename T, typename Basis, typename Problem, typename BilinearForm, typename RHSIntegral>
flens::DenseVector<flens::Array<T> > 
ThetaScheme<T, Basis, Problem, BilinearForm, RHSIntegral>::solve(T time, flens::DenseVector<flens::Array<T> > u_init, int level)
{
    op_LHSMatrix.setTime(time);
    op_RHSMatrix.setTime(time);
    op_RHSVector.setTime(time);
    flens::SparseGeMatrix<flens::CRS<T,flens::CRS_General> > lhsmatrix = problem.getStiffnessMatrix(op_LHSMatrix, level);
    flens::SparseGeMatrix<flens::CRS<T,flens::CRS_General> > rhsmatrix = problem.getStiffnessMatrix(op_RHSMatrix, level);
    flens::DenseVector<flens::Array<T> > rhsvector = problem.getRHS(op_RHSVector, level);
    flens::DenseVector<flens::Array<T> > rhsvector_total = rhsmatrix * u_init + rhsvector;
    flens::DenseVector<flens::Array<T> > u;
    cg(lhsmatrix, u, rhsvector_total);
    
    std::cout << "u(" << time << "): " << u << std::endl; 
    return u;
} 
  
}