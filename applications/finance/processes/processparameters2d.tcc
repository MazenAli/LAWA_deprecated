namespace lawa {

template <typename T>
ProcessParameters2D<T,BlackScholes2D>::ProcessParameters2D
(T _r, T _sigma1, T _sigma2, T _rho, T _u11, T _u12, T _u21, T _u22)
: r(_r), sigma1(_sigma1), sigma2(_sigma2), rho(_rho), u11(_u11), u12(_u12), u21(_u21), u22(_u22)
{

}

}   // namespace lawa
