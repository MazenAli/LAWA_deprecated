#ifndef LAWA_OPERATORS_HELMHOLTZOPERATOR1D_H
#define LAWA_OPERATORS_HELMHOLTZOPERATOR1D_H 1

#include <lawa/adaptive/index.h>
#include <lawa/integrals.h>
#include <lawa/enum.h>

namespace lawa {
    
/* HELMHOLTZ OPERATOR 
 *
 *    a(u,v) =  Integral(u_x * v_x) + c * Integral(u * v)
 *
 */
template <typename T, typename Basis>
class HelmholtzOperator1D{
    
    private:
        
        const Basis& basis;
        const T c;
        
        typedef typename Basis::BSplineType PrimalSpline;
        typedef typename Basis::WaveletType PrimalWavelet;
        
        PrimalSpline phi, d_phi;
        PrimalWavelet psi, d_psi;
        
        Integral<T, Gauss, PrimalSpline, PrimalSpline> integral_sfsf, dd_integral_sfsf;
        Integral<T, Gauss, PrimalSpline, PrimalWavelet> integral_sfw, dd_integral_sfw;
        Integral<T, Gauss, PrimalWavelet, PrimalSpline> integral_wsf, dd_integral_wsf;
        Integral<T, Gauss, PrimalWavelet, PrimalWavelet> integral_ww, dd_integral_ww;
            
    public:
        HelmholtzOperator1D(const Basis& _basis, const T _c);
        HelmholtzOperator1D(const HelmholtzOperator1D<T,Basis> &a);

        T getc() const;
        const Basis& getBasis() const;
    
        T
        operator()(XType xtype1, int j1, int k1, 
                   XType xtype2, int j2, int k2) const;

        T
        operator()(const Index1D &row_index, const Index1D &col_index) const;
    
        T
        operator()(T time, 
                   XType xtype1, int j1, int k1, 
                   XType xtype2, int j2, int k2) const
        {
            return operator()(xtype1, j1, k1, xtype2, j2, k2);
        }
};
    
    
} // namespace lawa

#include <lawa/operators/helmholtzoperator1d.tcc>

#endif // LAWA_OPERATORS_HELMHOLTZOPERATOR1D_H
