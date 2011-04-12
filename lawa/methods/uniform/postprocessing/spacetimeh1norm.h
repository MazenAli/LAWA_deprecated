#ifndef LAWA_METHODS_UNIFORM_POSTPROCESSING_SPACETIMEH1NORM_H
#define LAWA_METHODS_UNIFORM_POSTPROCESSING_SPACETIMEH1NORM_H 1
    
template <typename T, typename Basis>
class SpaceTimeH1Norm{
    
    private:
        
        const Basis& basis;
        
        typedef typename Basis::FirstBasisType::BSplineType PrimalSpline_t;
        typedef typename Basis::SecondBasisType::BSplineType PrimalSpline_x;
        typedef typename Basis::FirstBasisType::WaveletType PrimalWavelet_t;
        typedef typename Basis::SecondBasisType::WaveletType PrimalWavelet_x;
        
        PrimalSpline_t phi_t, d_phi_t;
        PrimalSpline_x phi_x, d_phi_x;
        PrimalWavelet_t psi_t, d_psi_t;
        PrimalWavelet_x psi_x, d_psi_x;
        
        lawa::Integral<T, lawa::Gauss, PrimalSpline_t, PrimalSpline_t> integral_sfsf_t,
                                                        dd_integral_sfsf_t;
        lawa::Integral<T, lawa::Gauss, PrimalSpline_x, PrimalSpline_x> integral_sfsf_x,
                                                        dd_integral_sfsf_x;
        lawa::Integral<T, lawa::Gauss, PrimalWavelet_t, PrimalWavelet_t>    integral_ww_t,
                                                             dd_integral_ww_t;
        lawa::Integral<T, lawa::Gauss, PrimalWavelet_x, PrimalWavelet_x>    integral_ww_x,
                                                             dd_integral_ww_x;
            
    public:
        SpaceTimeH1Norm(const Basis& _basis);
    
        T
        operator()(bool XisSpline, int j_t, int k_t, 
                   bool YisSpline, int j_x, int k_x) const;
};
    
    
#include "lawa/methods/uniform/postprocessing/spacetimeh1norm.tcc"

#endif // LAWA_METHODS_UNIFORM_POSTPROCESSING_SPACETIMEH1NORM_H