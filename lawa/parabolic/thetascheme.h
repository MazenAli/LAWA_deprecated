#ifndef LAWA_PARABOLIC_THETASCHEME_H
#define LAWA_PARABOLIC_THETASCHEME_H 1

#include <lawa/enum.h>

namespace lawa{
    
template<typename T, typename Basis, typename Problem, typename BilinearForm, typename RHSIntegral>
class ThetaScheme
{
    private:
        class Operator_LHSMatrix{
            private:
                ThetaScheme<T, Basis, Problem, BilinearForm, RHSIntegral>* scheme;
                const Basis& basis;
                const BilinearForm& a;
                T timestep;
                T time;
            
            public:                
                Operator_LHSMatrix(ThetaScheme<T, Basis, Problem, BilinearForm, RHSIntegral>* _scheme, const Basis& _basis, const BilinearForm& _a, T _timestep);
                
                T 
                operator()(XType xtype1, int j1, int k1,
                           XType xtype2, int j2, int k2) const;

                
                void setTime(T t){ time = t;}
            
        };
        
        class Operator_RHSMatrix{
            private:
                const ThetaScheme<T, Basis, Problem, BilinearForm, RHSIntegral>* scheme; 
                const Basis& basis;
                const BilinearForm& a;
                T timestep;
                T time;
            
            public:                
                Operator_RHSMatrix(const ThetaScheme<T, Basis, Problem, BilinearForm, RHSIntegral>* _scheme, const Basis& _basis, const BilinearForm& _a, T _timestep);
                
                T 
                operator()(XType xtype1, int j1, int k1,
                           XType xtype2, int j2, int k2) const;
                
                void setTime(T t){ time = t;}
            
            
        };
        
        class Operator_RHSVector{
            private:
                const ThetaScheme<T, Basis, Problem, BilinearForm, RHSIntegral>* scheme; 
                const Basis& basis;
                const RHSIntegral& rhs;
                T timestep;
                T time;
                
            public:                
                Operator_RHSVector(const ThetaScheme<T, Basis, Problem, BilinearForm, RHSIntegral>* _scheme, const Basis& _basis, const RHSIntegral& _rhs, T _timestep);
                
                T operator()(XType xtype, int j, int k) const;
                
                void setTime(T t){ time = t;}
                
            
        };
        
        friend class Operator_LHSMatrix;
        friend class Operator_RHSMatrix;
        friend class Operator_RHSVector;
        
        T theta;
        const Basis& basis;
        Problem problem;
        
        typedef typename Basis::BSplineType PrimalSpline;
        typedef typename Basis::WaveletType PrimalWavelet;
        PrimalSpline phi;
        PrimalWavelet psi;
        Integral<T, Gauss, PrimalSpline, PrimalSpline> integral_sfsf;
        Integral<T, Gauss, PrimalSpline, PrimalWavelet> integral_sfw;
        Integral<T, Gauss, PrimalWavelet, PrimalSpline> integral_wsf;
        Integral<T, Gauss, PrimalWavelet, PrimalWavelet> integral_ww;
        
        Operator_LHSMatrix op_LHSMatrix;
        Operator_RHSMatrix op_RHSMatrix;
        Operator_RHSVector op_RHSVector;

        
    public:        
        
        ThetaScheme(const T _theta, const Basis& _basis, const BilinearForm& _a, const RHSIntegral& _rhs, T _timestep);
    
        flens::DenseVector<flens::Array<T> > 
        solve(T time, flens::DenseVector<flens::Array<T> > u_init, int level);
        
        // Adaptive Erweiterung: Timestep in jedem Lösungsschritt neu setzen,
        //flens::DenseVector<flens::Array<T> > 
        //solve(T time, flens::DenseVector<flens::Array<T> > u_init, int level, T timestep);

};
      
} // namespace lawa

#include <lawa/parabolic/thetascheme.tcc>

#endif // LAWA_PARABOLIC_THETASCHEME_H