#ifndef LAWA_PARABOLIC_PROBLEM1D_H
#define LAWA_PARABOLIC_PROBLEM1D_H 1

#include <extensions/extensions.h>

namespace lawa{    

template<typename T, typename Basis>
class Problem1D
{
    private:
        const Basis& basis;
   
    public: 
        Problem1D(const Basis& _basis);
    
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

#include <lawa/parabolic/problem1d.tcc>

#endif // LAWA_PARABOLIC_PROBLEM1D_H