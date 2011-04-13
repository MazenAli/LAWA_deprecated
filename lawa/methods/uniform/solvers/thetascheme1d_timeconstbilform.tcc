namespace lawa{

// THETASCHEME
template<typename T, typename Basis, typename BilinearForm, typename RHSIntegral>
ThetaScheme1D_TimeConstBilForm<T, Basis, BilinearForm, RHSIntegral>::
ThetaScheme1D_TimeConstBilForm(const T _theta, const Basis& _basis, const BilinearForm& _a, RHSIntegral& _rhs)
    : theta(_theta), basis(_basis), problem(basis), phi(basis.mra), psi(basis),
      integral_sfsf(phi, phi), integral_sfw(phi, psi), integral_wsf(psi, phi), integral_ww(psi, psi),
      op_LHSMatrix(this, _a), op_RHSMatrix(this, _a), op_RHSVector(this, _rhs), prec(op_LHSMatrix),
      currentLevel(-1), P(problem.getPreconditioner(prec, 2))
{   
}
 
template<typename T, typename Basis, typename BilinearForm, typename RHSIntegral>
flens::DenseVector<flens::Array<T> > 
ThetaScheme1D_TimeConstBilForm<T, Basis, BilinearForm, RHSIntegral>::
solve(T time_old, T time_new, flens::DenseVector<flens::Array<T> > u_init, int level)
{   
    op_RHSVector.setTimes(time_old, time_new);
    if(level != currentLevel){
        op_LHSMatrix.setTimes(time_old, time_new);
        op_RHSMatrix.setTimes(time_old, time_new);
        lhsmatrix = problem.getStiffnessMatrix(op_LHSMatrix, level);
        rhsmatrix = problem.getStiffnessMatrix(op_RHSMatrix, level);
        P = problem.getPreconditioner(prec, level);
        currentLevel = level;
    }
    flens::DenseVector<flens::Array<T> > rhsvector = problem.getRHS(op_RHSVector, level);
    flens::DenseVector<flens::Array<T> > rhs = rhsmatrix * u_init + rhsvector;
    flens::DenseVector<flens::Array<T> > u(basis.mra.rangeI(level));
    pcg(P,lhsmatrix, u, rhs);
    //std::cout << cg(lhsmatrix, u, rhs) << "cg iterations" << std::endl;
    //std::cout << pcg(P, lhsmatrix, u, rhs) << "pcg iterations" << std::endl;
    
    //std::cout << "u(" << time_new << "): " << u << std::endl; 
    return u;
}

template<typename T, typename Basis, typename BilinearForm, typename RHSIntegral>
flens::DenseVector<flens::Array<T> > 
ThetaScheme1D_TimeConstBilForm<T, Basis, BilinearForm, RHSIntegral>::
solve(T time_old, T time_new, flens::DenseVector<flens::Array<T> > u_init, 
      flens::DenseVector<flens::Array<T> > f, int level)
{
     if(level != currentLevel){
         op_LHSMatrix.setTimes(time_old, time_new);
         op_RHSMatrix.setTimes(time_old, time_new);
         lhsmatrix = problem.getStiffnessMatrix(op_LHSMatrix, level);
         rhsmatrix = problem.getStiffnessMatrix(op_RHSMatrix, level);
         P = problem.getPreconditioner(prec, level);
         currentLevel = level;
     }
     flens::DenseVector<flens::Array<T> > rhs = rhsmatrix * u_init + f;
     flens::DenseVector<flens::Array<T> > u(u_init);
     pcg(P, lhsmatrix, u, rhs);
     //std::cout << cg(lhsmatrix, u, rhs) << "cg iterations" << std::endl;
     //std::cout << pcg(P, lhsmatrix, u, rhs) << "pcg iterations" << std::endl;
     //std::cout << "u(" << time_new << "): " << u << std::endl; 
     return u;     
}


template<typename T, typename Basis, typename BilinearForm, typename RHSIntegral>
void
ThetaScheme1D_TimeConstBilForm<T, Basis, BilinearForm, RHSIntegral>::
setRHS(RHSIntegral& _rhs)
{
    op_RHSVector.setRHS(_rhs);
}

template<typename T, typename Basis, typename BilinearForm, typename RHSIntegral>
flens::SparseGeMatrix<flens::CRS<T,flens::CRS_General> > 
ThetaScheme1D_TimeConstBilForm<T, Basis, BilinearForm, RHSIntegral>::
getLHSMatrix(int level)
{   
    return lhsmatrix;
}


/*======================================================================================*/    
// OPERATOR_LHSMATRIX
template<typename T, typename Basis, typename BilinearForm, typename RHSIntegral>
ThetaScheme1D_TimeConstBilForm<T, Basis, BilinearForm, RHSIntegral>::Operator_LHSMatrix::
Operator_LHSMatrix(ThetaScheme1D_TimeConstBilForm<T, Basis, BilinearForm, RHSIntegral>* _scheme, 
                   const BilinearForm& _a)
    : a(_a)
{   
    scheme = _scheme;
}

template<typename T, typename Basis, typename BilinearForm, typename RHSIntegral>
T 
ThetaScheme1D_TimeConstBilForm<T, Basis, BilinearForm, RHSIntegral>::Operator_LHSMatrix::
operator()(XType xtype1, int j1, int k1,
           XType xtype2, int j2, int k2) const
{
    // (M + deltaT * theta * A_k+1)
    return scheme->integral(j1, k1, xtype1, 0, j2, k2, xtype2, 0) 
            + (time_new - time_old) * scheme->theta * a(xtype1, j1, k1, xtype2, j2, k2);
}


// OPERATOR_RHSMATRIX
template<typename T, typename Basis, typename BilinearForm, typename RHSIntegral>
ThetaScheme1D_TimeConstBilForm<T, Basis, BilinearForm, RHSIntegral>::Operator_RHSMatrix::
Operator_RHSMatrix(const ThetaScheme1D_TimeConstBilForm<T, Basis, BilinearForm, RHSIntegral>* _scheme, 
                   const BilinearForm& _a)
    : a(_a)
{
     scheme = _scheme;
}


template<typename T, typename Basis, typename BilinearForm, typename RHSIntegral>
T 
ThetaScheme1D_TimeConstBilForm<T, Basis, BilinearForm, RHSIntegral>::Operator_RHSMatrix::
operator()(XType xtype1, int j1, int k1,
           XType xtype2, int j2, int k2) const
{
   // (M - deltaT * (1-theta) * A_k)
   return scheme->integral(j1, k1, xtype1, 0, j2, k2, xtype2, 0)
        - (time_new - time_old) * (1. - scheme->theta) * a(xtype1,j1,k1, xtype2, j2,k2);
}

// OPERATOR_RHSVECTOR
template<typename T, typename Basis, typename BilinearForm, typename RHSIntegral>
ThetaScheme1D_TimeConstBilForm<T, Basis, BilinearForm, RHSIntegral>::Operator_RHSVector::
Operator_RHSVector(const ThetaScheme1D_TimeConstBilForm<T, Basis, BilinearForm, RHSIntegral>* _scheme, 
                   RHSIntegral& _rhs)
    : rhs(_rhs)
{
     scheme = _scheme;
}

template<typename T, typename Basis, typename BilinearForm, typename RHSIntegral>
T 
ThetaScheme1D_TimeConstBilForm<T, Basis, BilinearForm, RHSIntegral>::Operator_RHSVector::
operator()(XType xtype, int j, int k) const
{   
    // deltaT * (theta*f_k+1 - (1-theta)*f_k)
    return (time_new - time_old)*(scheme->theta * rhs(time_new, xtype, j, k) 
                    + (1. - scheme->theta)*rhs(time_old, xtype, j, k));
} 
  
}