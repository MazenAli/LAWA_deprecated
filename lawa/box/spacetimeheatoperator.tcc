namespace lawa{
    
template <typename T, typename Basis>    
SpaceTimeHeatOperator<T, Basis>::SpaceTimeHeatOperator(const Basis& _basis, const T _c)
    : basis(_basis), c(_c), 
      phi_t(basis.first.mra), d_phi_t(basis.first.mra, 1),
      phi_x(basis.second.mra), d_phi_x(basis.second.mra, 1),
      psi_t(basis.first), d_psi_t(basis.first, 1),
      psi_x(basis.second), d_psi_x(basis.second, 1),
      integral_sfsf_t(phi_t, phi_t), d_integral_sfsf_t(d_phi_t, phi_t),
      integral_sfsf_x(phi_x, phi_x), dd_integral_sfsf_x(d_phi_x, d_phi_x),
      integral_sfw_t(phi_t, psi_t), d_integral_sfw_t(d_phi_t, psi_t),
      integral_sfw_x(phi_x, psi_x), dd_integral_sfw_x(d_phi_x, d_psi_x),
      integral_wsf_t(psi_t, phi_t), d_integral_wsf_t(d_psi_t, phi_t),
      integral_wsf_x(psi_x, phi_x), dd_integral_wsf_x(d_psi_x, d_phi_x),
      integral_ww_t(psi_t, psi_t), d_integral_ww_t(d_psi_t, psi_t),
      integral_ww_x(psi_x, psi_x), dd_integral_ww_x(d_psi_x, d_psi_x)
{
}
    
template <typename T, typename Basis>      
T
SpaceTimeHeatOperator<T, Basis>::operator()(bool FirstXisSpline, bool FirstYisSpline,
                                        int j1_t, int k1_t, int j1_x, int k1_x,
                                        bool SecondXisSpline, bool SecondYisSpline,
                                        int j2_t, int k2_t, int j2_x, int k2_x) const
{
    T val_t = 0;
    T d_val_t = 0;
    T val_x = 0;
    T dd_val_x = 0;
    
    if(FirstXisSpline){
         if(SecondXisSpline){
             val_t = integral_sfsf_t(j1_t, k1_t, j2_t, k2_t);
             d_val_t = d_integral_sfsf_t(j1_t, k1_t, j2_t, k2_t);
         }
         else{
             val_t = integral_sfw_t(j1_t, k1_t, j2_t, k2_t);
             d_val_t = d_integral_sfw_t(j1_t, k1_t, j2_t, k2_t);
         }
    }
    else{
         if(SecondXisSpline){
             val_t = integral_wsf_t(j1_t, k1_t, j2_t, k2_t);
             d_val_t = d_integral_wsf_t(j1_t, k1_t, j2_t, k2_t);
         }
         else{
             val_t = integral_ww_t(j1_t, k1_t, j2_t, k2_t);
             d_val_t = d_integral_ww_t(j1_t, k1_t, j2_t, k2_t);
         }
    }
    
    if(FirstYisSpline){
         if(SecondYisSpline){
             val_x = integral_sfsf_x(j1_x, k1_x, j2_x, k2_x);
             dd_val_x = dd_integral_sfsf_x(j1_x, k1_x, j2_x, k2_x);
         }
         else{
             val_x = integral_sfw_x(j1_x, k1_x, j2_x, k2_x);
             dd_val_x = dd_integral_sfw_x(j1_x, k1_x, j2_x, k2_x);
         }
    }
    else{
         if(SecondYisSpline){
             val_x = integral_wsf_x(j1_x, k1_x, j2_x, k2_x);
             dd_val_x = dd_integral_wsf_x(j1_x, k1_x, j2_x, k2_x);
         }
         else{
             val_x = integral_ww_x(j1_x, k1_x, j2_x, k2_x);
             dd_val_x = dd_integral_ww_x(j1_x, k1_x, j2_x, k2_x);
         }
    }
    
    
    return d_val_t * val_x  + c * val_t * dd_val_x;
}
    
} // namespace lawa