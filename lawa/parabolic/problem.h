#ifndef LAWA_PARABOLIC_PROBLEM_H
#define LAWA_PARABOLIC_PROBLEM_H 1

#include <extensions/extensions.h>

namespace lawa{    

template<typename T, typename Basis>
class Problem
{
    private:
        const Basis& basis;
   
    public: 
        Problem(const Basis& _basis);
    
        template <typename BilinearForm>
        flens::SparseGeMatrix<flens::CRS<T,flens::CRS_General> >
        getStiffnessMatrix(BilinearForm& a, int J, T tol = 10e-15);
    
        template <typename RHSIntegral>
        flens::DenseVector<flens::Array<T> >
        getRHS(RHSIntegral& rhs, int J);
        
        template <typename Preconditioner>
        flens::DiagonalMatrix<T>    
        getPreconditioner(Preconditioner& P, int J);
};

} // namespace lawa

#include <lawa/parabolic/problem.tcc>

#endif // LAWA_PARABOLIC_PROBLEM_H