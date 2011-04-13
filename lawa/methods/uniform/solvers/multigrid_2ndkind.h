#ifndef LAWA_METHODS_UNIFORM_SOLVERS_MULTIGRID_2NDKIND_H
#define LAWA_METHODS_UNIFORM_SOLVERS_MULTIGRID_2NDKIND_H 1

#include <lawa/methods/uniform/datastructures/datastructures.h>
#include <lawa/methods/uniform/solvers/thetascheme1d.h>
#include <lawa/righthandsides/homogeneousrhs.h>

namespace lawa{
    
template<typename T, typename PrimalBasis, typename DualBasis, typename BilinearForm, typename RHSIntegral>
class MultiGrid_2ndKind{
    
    public: 
        
        MultiGrid_2ndKind(PrimalBasis& _b, DualBasis& _b_, BilinearForm& a, RHSIntegral& rhs, 
                          T theta, T deltaT, int timesteps, int minLevel);
        
        flens::DenseVector<flens::Array<T> > 
        run_MG_2ndKind(flens::DenseVector<flens::Array<T> >& u0, int maxLevel);
    
    private:
        
        PrimalBasis& b;
        DualBasis& b_;
        typedef ThetaScheme1D<T, PrimalBasis, BilinearForm, RHSIntegral> FullThetaScheme;
        typedef TimeStepping<T, FullThetaScheme> FullTimeStepMethod;    
        typedef ThetaScheme1D<T, PrimalBasis, BilinearForm, HomogeneousRHS<T> > HomThetaScheme;
        typedef TimeStepping<T, HomThetaScheme> HomTimeStepMethod;
        typedef FixedPointSolver<T, HomTimeStepMethod> ThetaFPSolver;
        
        FullThetaScheme full_theta;
        FullTimeStepMethod full_ts;
        HomogeneousRHS<T> hom_rhs;
        HomThetaScheme hom_theta;
        HomTimeStepMethod hom_ts;
        ThetaFPSolver fp;


        class MG_2ndKind_Smoother{
            private:
                MultiGrid_2ndKind<T, PrimalBasis, DualBasis, BilinearForm, RHSIntegral>* mg_ptr;
        
            public: 
                MG_2ndKind_Smoother(MultiGrid_2ndKind<T, PrimalBasis, DualBasis, BilinearForm, RHSIntegral>* ref);
                
                flens::DenseVector<flens::Array<T> > 
                solve(flens::DenseVector<flens::Array<T> > u, flens::DenseVector<flens::Array<T> > f);  
                
                void setLevel(int level){ mg_ptr->hom_ts.setLevel(level);}
                
                flens::DenseVector<flens::Array<T> >
                getResiduum(flens::DenseVector<flens::Array<T> >& u,flens::DenseVector<flens::Array<T> >& f)
                { return mg_ptr->hom_ts.getResiduum(u) - f;}
        };
    
        class MG_2ndKind_Solver{
            private:
                MultiGrid_2ndKind<T, PrimalBasis, DualBasis, BilinearForm, RHSIntegral>* mg_ptr;
        
            public: 
                MG_2ndKind_Solver(MultiGrid_2ndKind<T, PrimalBasis, DualBasis, BilinearForm, RHSIntegral>* ref);
                
                flens::DenseVector<flens::Array<T> > 
                solve(flens::DenseVector<flens::Array<T> > u, flens::DenseVector<flens::Array<T> > f);  
                
                void setLevel(int level){ mg_ptr->fp.setLevel(level);}
        };
        
        friend class MG_2ndKind_Smoother;
        friend class MG_2ndKind_Solver;
        
        typedef MultigridSolver<T, PrimalBasis, DualBasis, MG_2ndKind_Smoother, MG_2ndKind_Solver> Multigrid;
        MG_2ndKind_Smoother mg_smoother;
        MG_2ndKind_Solver mg_solver;
        Multigrid mg;

};



} // namespace lawa

#include <lawa/methods/uniform/solvers/multigrid_2ndkind.tcc>

#endif // LAWA_METHODS_UNIFORM_SOLVERS_MULTIGRID_2NDKIND_H