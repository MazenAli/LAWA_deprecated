namespace lawa {



template <typename SingularKernel, typename First, typename Second>
typename RestrictTo<!IsDual<First>::value and !IsDual<Second>::value, typename First::T>::Type
_integrate(const SingularIntegral<SingularKernel,First,Second> &singularintegral)
{
    typedef typename First::T T;
    // note: called phi_... but can also be a wavelet.
    Support<T> supp_x = singularintegral.first.generator(singularintegral.e1).support(singularintegral.j1,singularintegral.k1);
    Support<T> supp_y = singularintegral.second.generator(singularintegral.e2).support(singularintegral.j2,singularintegral.k2);

    DenseVector<Array<T> > SingularPoints_x =
             singularintegral.first.generator(singularintegral.e1).singularSupport(singularintegral.j1,singularintegral.k1);
    DenseVector<Array<T> > SingularPoints_y =
            singularintegral.second.generator(singularintegral.e2).singularSupport(singularintegral.j2,singularintegral.k2);

    std::vector<T> AllSingularPoints_x, AllSingularPoints_y;
    for (int i=SingularPoints_x.firstIndex(); i<=SingularPoints_x.lastIndex(); ++i) {
        T x_i = SingularPoints_x(i);
        AllSingularPoints_x.push_back(x_i);
        if (supp_y.l1 <= x_i && x_i <= supp_y.l2) AllSingularPoints_y.push_back(x_i);
    }
    for (int i=SingularPoints_y.firstIndex(); i<=SingularPoints_y.lastIndex(); ++i) {
        T y_i = SingularPoints_y(i);
        AllSingularPoints_y.push_back(y_i);
        if (supp_x.l1 <= y_i && y_i <= supp_x.l2) AllSingularPoints_x.push_back(y_i);
    }

    sort(AllSingularPoints_x.begin(),  AllSingularPoints_x.end() );
    typename std::vector<T>::iterator it_x;
    it_x = unique(AllSingularPoints_x.begin(),  AllSingularPoints_x.end() );
    AllSingularPoints_x.resize(it_x - AllSingularPoints_x.begin());

    typename std::vector<T>::iterator it_y;
    sort(AllSingularPoints_y.begin(),  AllSingularPoints_y.end() );
    it_y = unique(AllSingularPoints_y.begin(),  AllSingularPoints_y.end() );
    AllSingularPoints_y.resize(it_y - AllSingularPoints_y.begin());


    typename std::vector<T>::const_iterator it_x1, it_x2, it_y1, it_y2;
    it_x1=AllSingularPoints_x.begin();
    it_x2=AllSingularPoints_x.begin(); ++it_x2;
    it_y1=AllSingularPoints_y.begin();
    it_y2=AllSingularPoints_y.begin(); ++it_y2;

    long double ret = 0.L;
    while (it_x2!=AllSingularPoints_x.end()) {
        long double a1 = (*it_x1), b1 = (*it_x2);
        while (it_y2!=AllSingularPoints_y.end()) {
            long double a2 = (*it_y1), b2 = (*it_y2);
            ret += singularintegral.singularquadrature(a1, b1, a2, b2, (long double)1e-10);
            ++it_y1; ++it_y2;
        }
        ++it_x1; ++it_x2;
    }
    return (T)ret;
}


template <typename SingularKernel, typename First, typename Second>
SingularIntegral<SingularKernel,First,Second>::SingularIntegral(const SingularKernel &_singularkernel,
                                                                const First &_first,
                                                                const Second &_second)
    : singularkernel(_singularkernel), first(_first), second(_second), singularquadrature(*this)
{
}

template <typename SingularKernel, typename First, typename Second>
typename First::T
SingularIntegral<SingularKernel,First,Second>::p1(T x) const
{
    return first.generator(e1).operator()(x,j1,k1,deriv1);
}

template <typename SingularKernel, typename First, typename Second>
typename First::T
SingularIntegral<SingularKernel,First,Second>::p2(T x) const
{
    return second.generator(e2).operator()(x,j2,k2,deriv2);
}

template <typename SingularKernel, typename First, typename Second>
typename First::T
SingularIntegral<SingularKernel,First,Second>::kernel(T x) const
{
    return singularkernel.operator()(x);
}

template <typename SingularKernel, typename First, typename Second>
typename First::T
SingularIntegral<SingularKernel,First,Second>::operator()
                                               (int _j1, long _k1, XType _e1, int _deriv1,
                                                int _j2, long _k2, XType _e2, int _deriv2) const
{
    j1 = _j1; k1 = _k1; e1 = _e1; deriv1 = _deriv1;
    j2 = _j2; k2 = _k2; e2 = _e2; deriv2 = _deriv2;
    return _integrate(*this);
}

}   // namespace lawa

