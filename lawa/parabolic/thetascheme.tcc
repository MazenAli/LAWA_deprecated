namespace lawa{

    
// OPERATOR_LHSMATRIX
template<typename T, typename Basis, typename Problem, typename BilinearForm, typename RHSIntegral>
ThetaScheme<T, Basis, Problem, BilinearForm, RHSIntegral>::Operator_LHSMatrix::
Operator_LHSMatrix(const Basis& _basis, const BilinearForm& _a, T _timestep)
    : basis(_basis), a(_a), timestep(_timestep), time(0), phi(basis.mra), psi(basis),
      integral_sfsf(phi, phi), integral_sfw(phi, psi), integral_wsf(psi, phi), integral_ww(psi, psi)
{
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
            return integral_sfsf(j1,k1,j2,k2) + timestep*theta*a(time, xtype1,j1,k1, xtype2, j2,k2);
        }
        else{
            return integral_sfw(j1,k1,j2,k2) + timestep*theta*a(time, xtype1,j1,k1, xtype2, j2,k2);            
        }
    }
    else{
        if(xtype2 == XBSpline){
            return integral_wsf(j1,k1,j2,k2) + timestep*theta*a(time, xtype1,j1,k1, xtype2, j2,k2);
        }
        else{
            return integral_ww(j1,k1,j2,k2) + timestep*theta*a(time, xtype1,j1,k1, xtype2, j2,k2);
        }
    }
    
}

// OPERATOR_RHSMATRIX
template<typename T, typename Basis, typename Problem, typename BilinearForm, typename RHSIntegral>
ThetaScheme<T, Basis, Problem, BilinearForm, RHSIntegral>::Operator_RHSMatrix::
Operator_RHSMatrix(const Basis& _basis, const BilinearForm& _a, T _timestep)
    : basis(_basis), a(_a), timestep(_timestep), time(0), phi(basis.mra), psi(basis),
      integral_sfsf(phi, phi), integral_sfw(phi, psi), integral_wsf(psi, phi), integral_ww(psi, psi)
{
}

// OPERATOR_RHSVECTOR
template<typename T, typename Basis, typename Problem, typename BilinearForm, typename RHSIntegral>
ThetaScheme<T, Basis, Problem, BilinearForm, RHSIntegral>::Operator_RHSVector::
Operator_RHSVector(const Basis& _basis, const RHSIntegral& _rhs, T _timestep)
    : basis(_basis), rhs(_rhs), timestep(_timestep), time(0) 
{
}


// THETASCHEME
template<typename T, typename Basis, typename Problem, typename BilinearForm, typename RHSIntegral>
ThetaScheme<T, Basis, Problem, BilinearForm, RHSIntegral>::
ThetaScheme(const T _theta, const Basis& _basis, const BilinearForm& _a, const RHSIntegral& _rhs, T _timestep)
    : basis(_basis), problem(basis), op_LHSMatrix(basis, _a, _timestep),
      op_RHSMatrix(basis,_a, _timestep), op_RHSVector(basis, _rhs, _timestep)
{
}
 
template<typename T, typename Basis, typename Problem, typename BilinearForm, typename RHSIntegral>
flens::DenseVector<flens::Array<T> > 
ThetaScheme<T, Basis, Problem, BilinearForm, RHSIntegral>::solve(T time, flens::DenseVector<flens::Array<T> > u_init, int level)
{
    flens::DenseVector<flens::Array<T> > u;
    std::cout << op_LHSMatrix(XBSpline, basis.j0, 1, XBSpline, basis.j0, 1) << std::endl;
    return u;
} 
  
}