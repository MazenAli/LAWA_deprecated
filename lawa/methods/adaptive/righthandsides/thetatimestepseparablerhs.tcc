namespace lawa {

template <typename T, typename Index, typename SpatialRHS>
ThetaTimeStepSeparableRHS<T,Index,SpatialRHS>::ThetaTimeStepSeparableRHS
(Function<T> &_fct_t, SpatialRHS &_F_x)
: fct_t(_fct_t), F_x(_F_x),
  discrete_timepoint(0.), theta(1.), timestep(0.1)
{

}

template <typename T, typename Index, typename SpatialRHS>
void
ThetaTimeStepSeparableRHS<T,Index,SpatialRHS>::setThetaTimeStepParameters
(T _theta,T _timestep, T _discrete_timepoint, Coefficients<Lexicographical,T,Index> *_propagated_u_k)
{
    theta = _theta;
    timestep = _timestep;
    discrete_timepoint = _discrete_timepoint;
    propagated_u_k = _propagated_u_k;
    //std::cerr << "propagated_u_k = " << *propagated_u_k << std::endl;
}

template <typename T, typename Index, typename SpatialRHS>
T
ThetaTimeStepSeparableRHS<T,Index,SpatialRHS>::operator()(T t, const Index &index)
{
    return fct_t(t) * F_x(index);
}

template <typename T, typename Index, typename SpatialRHS>
T
ThetaTimeStepSeparableRHS<T,Index,SpatialRHS>::operator()(const Index &index)
{
    T spatial_val = F_x(index);
    coeff_it it = (*propagated_u_k).find(index);

    T return_val =  timestep*(   fct_t(discrete_timepoint)         *   theta  *spatial_val );
    if (theta != 1.)  {
        return_val  += timestep*( fct_t(discrete_timepoint-timestep)*(1.-theta)*spatial_val );
    }
    if (it!=(*propagated_u_k).end()) return_val += (*it).second;
    return return_val;
}

}   // namespace lawa
