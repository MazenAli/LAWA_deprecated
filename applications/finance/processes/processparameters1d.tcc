namespace lawa {

template <typename T>
ProcessParameters1D<T,CGMY>::ProcessParameters1D(T _r, T _k_C, T _k_G, T _k_M, T _k_Y)
: r(_r), k_C(_k_C), k_G(_k_G), k_M(_k_M), k_Y(_k_Y)
{

}

template <typename T>
ProcessParameters1D<T,BlackScholes>::ProcessParameters1D(T _r, T _sigma)
: r(_r), sigma(_sigma)
{

}

}   // namespace lawa
