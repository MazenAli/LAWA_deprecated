namespace lawa {

template<typename T>
DiagonalScalingPreconditioner2D<T>::DiagonalScalingPreconditioner2D(int _s_x, int _s_y)
    : s_x(_s_x), s_y(_s_y)
{
    assert(s_x >= 0);
    assert(s_y >= 0);
}

template<typename T>
T
DiagonalScalingPreconditioner2D<T>::operator()(bool XisSpline, int j_x, int k_x,
                                               bool YisSpline, int j_y, int k_y) const
{
    return pow2i<T>(- s_x*j_x - s_y*j_y);
}



} // namespace lawa
