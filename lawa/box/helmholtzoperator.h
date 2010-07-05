#ifndef LAWA_BOX_HELMHOLTZOPERATOR_H
#define LAWA_BOX_HELMHOLTZOPERATOR_H 1


#include <lawa/integrals.h>

namespace lawa {
    
template <typename T, typename Basis>
class HelmholtzOperator{
    
    private:
        
        Basis basis;
        T c;
        
        typedef typename Basis::FirstBasisType::BSplineType PrimalSpline_x;
        typedef typename Basis::SecondBasisType::BSplineType PrimalSpline_y;
        typedef typename Basis::FirstBasisType::WaveletType PrimalWavelet_x;
        typedef typename Basis::SecondBasisType::WaveletType PrimalWavelet_y;
        
        PrimalSpline_x phi_x, d_phi_x;
        PrimalSpline_y phi_y, d_phi_y;
        PrimalWavelet_x psi_x, d_psi_x;
        PrimalWavelet_y psi_y, d_psi_y;
        
        Integral<T, Gauss, PrimalSpline_x, PrimalSpline_x> integral_sfsf_x,
                                                        dd_integral_sfsf_x;
        Integral<T, Gauss, PrimalSpline_y, PrimalSpline_y> integral_sfsf_y,
                                                        dd_integral_sfsf_y;
        Integral<T, Gauss, PrimalSpline_x, PrimalWavelet_x>    integral_sfw_x,
                                                            dd_integral_sfw_x;
        Integral<T, Gauss, PrimalSpline_y, PrimalWavelet_y>    integral_sfw_y,
                                                            dd_integral_sfw_y;                                      
        Integral<T, Gauss, PrimalWavelet_x, PrimalSpline_x>    integral_wsf_x,
                                                            dd_integral_wsf_x;
        Integral<T, Gauss, PrimalWavelet_y, PrimalSpline_y>    integral_wsf_y,
                                                            dd_integral_wsf_y;
        Integral<T, Gauss, PrimalWavelet_x, PrimalWavelet_x>    integral_ww_x,
                                                             dd_integral_ww_x;
        Integral<T, Gauss, PrimalWavelet_y, PrimalWavelet_y>    integral_ww_y,
                                                             dd_integral_ww_y;
            
    public:
        HelmholtzOperator(Basis _basis, T _c);
    
        T
        operator()(bool FirstXisSpline, bool FirstYisSpline,
                   bool SecondXisSpline, bool SecondYisSpline,
                   int j1_x, int k1_x, int j2_x, int k2_x,
                   int j1_y, int k1_y, int j2_y, int k2_y) const;
    
};
    
    
} // namespace lawa

#include <lawa/box/helmholtzoperator.tcc>

#endif // LAWA_BOX_HELMHOLTZOPERATOR_H