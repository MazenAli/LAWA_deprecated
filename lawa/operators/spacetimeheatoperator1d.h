#ifndef LAWA_OPERATORS_SPACETIMEHEATOPERATOR1D_H
#define LAWA_OPERATORS_SPACETIMEHEATOPERATOR1D_H 1


#include <lawa/integrals.h>

namespace lawa {
    
template <typename T, typename Basis>
class SpaceTimeHeatOperator1D{
    
    private:
        
        const Basis& basis;
        const T c;
        
        typedef typename Basis::FirstBasisType::BSplineType PrimalSpline_t;
        typedef typename Basis::SecondBasisType::BSplineType PrimalSpline_x;
        typedef typename Basis::FirstBasisType::WaveletType PrimalWavelet_t;
        typedef typename Basis::SecondBasisType::WaveletType PrimalWavelet_x;
        
        PrimalSpline_t phi_t, d_phi_t;
        PrimalSpline_x phi_x, d_phi_x;
        PrimalWavelet_t psi_t, d_psi_t;
        PrimalWavelet_x psi_x, d_psi_x;
        
        Integral<T, Gauss, PrimalSpline_t, PrimalSpline_t> integral_sfsf_t,
                                                        d_integral_sfsf_t;
        Integral<T, Gauss, PrimalSpline_x, PrimalSpline_x> integral_sfsf_x,
                                                        dd_integral_sfsf_x;
        Integral<T, Gauss, PrimalSpline_t, PrimalWavelet_t>    integral_sfw_t,
                                                            d_integral_sfw_t;
        Integral<T, Gauss, PrimalSpline_x, PrimalWavelet_x>    integral_sfw_x,
                                                            dd_integral_sfw_x;                                      
        Integral<T, Gauss, PrimalWavelet_t, PrimalSpline_t>    integral_wsf_t,
                                                            d_integral_wsf_t;
        Integral<T, Gauss, PrimalWavelet_x, PrimalSpline_x>    integral_wsf_x,
                                                            dd_integral_wsf_x;
        Integral<T, Gauss, PrimalWavelet_t, PrimalWavelet_t>    integral_ww_t,
                                                             d_integral_ww_t;
        Integral<T, Gauss, PrimalWavelet_x, PrimalWavelet_x>    integral_ww_x,
                                                             dd_integral_ww_x;
            
    public:
        SpaceTimeHeatOperator1D(const Basis& _basis, const T _c);
        
        T getc() const;
        const Basis& getBasis() const;
    
        T                                                           // returns a(v,u)
        operator()(XType row_xtype_t, int j1_t, int k1_t, 
                   XType row_xtype_x, int j1_x, int k1_x,
                   XType col_xtype_t,int j2_t, int k2_t, 
                   XType col_xtype_x,int j2_x, int k2_x) const;
        T
        operator()(XType xtype_t, int j_t, int k_t,                 // returns a(u,u)
                   XType xtype_x, int j_x, int k_x) const;
        
        T
        operator()(const Index2D &row_index, const Index2D &col_index) const;
};
    
    
} // namespace lawa

#include <lawa/operators/spacetimeheatoperator1d.tcc>

#endif // LAWA_OPERATORS_SPACETIMEHEATOPERATOR1D_H