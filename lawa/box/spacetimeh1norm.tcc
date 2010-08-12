#include <iostream>

template<typename T, typename Basis>
SpaceTimeH1Norm<T,Basis>::SpaceTimeH1Norm(const Basis& _basis)
    : basis(_basis),
      phi_t(basis.first.mra), d_phi_t(basis.first.mra, 1),
      phi_x(basis.second.mra), d_phi_x(basis.second.mra, 1),
      psi_t(basis.first), d_psi_t(basis.first, 1),
      psi_x(basis.second), d_psi_x(basis.second, 1),
      integral_sfsf_t(phi_t, phi_t), dd_integral_sfsf_t(d_phi_t, d_phi_t),
      integral_sfsf_x(phi_x, phi_x), dd_integral_sfsf_x(d_phi_x, d_phi_x),
      integral_ww_t(psi_t, psi_t), dd_integral_ww_t(d_psi_t, d_psi_t),
      integral_ww_x(psi_x, psi_x), dd_integral_ww_x(d_psi_x, d_psi_x){
}

template<typename T, typename Basis>
T
SpaceTimeH1Norm<T,Basis>::operator()(bool XisSpline, int j_t, int k_t, 
                                       bool YisSpline, int j_x, int k_x) const
{
    T norm_t_L2 = 0;
    //T norm_t_H1 = 0;
    T norm_x_H1 = 0;
    
    if(XisSpline){
        norm_t_L2 = integral_sfsf_t(j_t, k_t, j_t, k_t);
        //norm_t_H1 = norm_t_L2 + dd_integral_sfsf_t(j_t, k_t, j_t, k_t);
        //std::cout << "Spline j=" << j_t <<", k=" << k_t << "  Norm_t = " << norm_t 
       // << " = " <<  integral_sfsf_t(j_t, k_t, j_t, k_t)<<  " + " << dd_integral_sfsf_t(j_t, k_t, j_t, k_t) << std::endl;
    }
    else{
        norm_t_L2 = integral_ww_t(j_t, k_t, j_t, k_t);
       // norm_t_H1 = norm_t_L2 + dd_integral_ww_t(j_t, k_t, j_t, k_t);
       // std::cout << "Wavelet j=" << j_t <<", k=" << k_t <<"  Norm_t = " << norm_t 
       // << " = " <<  integral_ww_t(j_t, k_t, j_t, k_t)<<  " + " << dd_integral_ww_t(j_t, k_t, j_t, k_t) << std::endl;
    }
  
    if(YisSpline){
        norm_x_H1 = integral_sfsf_x(j_x, k_x, j_x, k_x) + dd_integral_sfsf_x(j_x, k_x, j_x, k_x);
       // std::cout << "Spline j=" << j_x <<", k=" << k_x << "  Norm_xV = " << norm_xV 
       // << " = " <<  integral_sfsf_x(j_x, k_x, j_x, k_x)<<  " + " << dd_integral_sfsf_x(j_x, k_x, j_x, k_x)<< std::endl;
        //norm_xV_ = ????
    }
    else{
        norm_x_H1 = integral_ww_x(j_x, k_x, j_x, k_x) + dd_integral_ww_x(j_x, k_x, j_x, k_x);
       // std::cout << "Wavelet j=" << j_x <<", k=" << k_x <<"  Norm_xV = " << norm_xV
       // << " = " <<  integral_ww_x(j_x, k_x, j_x, k_x)<<  " + " << dd_integral_ww_x(j_x, k_x, j_x, k_x)<< std::endl;
    }
    
    
    return (norm_t_L2 * norm_x_H1 + lawa::pow2i<T>(2*j_t) * norm_t_L2 / norm_x_H1);                              
}