#ifndef LAWA_BOX_BOXPROBLEM_H
#define LAWA_BOX_BOXPROBLEM_H 1

namespace lawa{    

template<typename T, typename Basis, typename BilinearForm>
//template<typename T, typename Basis, typename BilinearForm, typename RHSIntegral, typename Preconditioner>
class BoxProblem
{
    private:
        Basis basis;
        BilinearForm a;
        //RHSIntegral rhs;
        //Preconditioner P;
   
    public: 
        BoxProblem(Basis _basis, BilinearForm _a);
        //BoxProblem(Basis _basis, BilinearForm _a, RHSIntegral _rhs, Preconditioner _P);
    
        flens::GeMatrix<flens::FullStorage<T, cxxblas::ColMajor> >
        getStiffnessMatrix(int J_x, int J_y, T tol = 10e-15);
    
        //flens::GeMatrix<flens::FullStorage<T, flens::ColMajor> >    
        //getPreconditioner(basis, P, int J_x, int J_y);
    
        //flens::DenseVector<flens::Array<T> >
        //getRHS(basis, rhs, int J_x, int J_y);

};

} // namespace lawa

#include <lawa/box/boxproblem.tcc>

#endif // LAWA_BOX_BOXPROBLEM_H