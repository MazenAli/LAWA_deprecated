namespace lawa {

template <OptionType1D OType, typename Basis1D>
typename Basis1D::T
_integrate(const PayoffInitialCondition1D<OType,Basis1D> &initcond)
{
    typedef typename Basis1D::T T;
    const typename Basis1D::BasisFunctionType &basisfunction = initcond.basis.generator(initcond.e1);

    DenseVector<Array<T> > singularPoints;

    int p = initcond.option.singularPoints.length();

    const DenseVector<Array<T> > & basisfunctionSingularPoints
                                    = basisfunction.singularSupport(initcond.j1,initcond.k1);
    int m = basisfunctionSingularPoints.length();
    singularPoints.engine().resize(m+p+1);

    std::merge(basisfunctionSingularPoints.engine().data(),
               basisfunctionSingularPoints.engine().data() + m,
               initcond.option.singularPoints.engine().data(),
               initcond.option.singularPoints.engine().data() + p,
               singularPoints.engine().data());

    singularPoints(m+p+1) = 0.;                         //singular point of weight
    std::inplace_merge(singularPoints.engine().data(),
                       singularPoints.engine().data() + m + p,
                       singularPoints.engine().data() + m + p + 1);

    T ret = 0.0;
    for (int i=singularPoints.firstIndex(); i<singularPoints.lastIndex(); ++i) {
        ret += initcond.quadrature(singularPoints(i),singularPoints(i+1));
    }
    return ret;
}

template <OptionType1D OType, typename Basis1D>
typename Basis1D::T
_integrand(const PayoffInitialCondition1D<OType,Basis1D> &initcond,
                            typename Basis1D::T x)
{
    const typename Basis1D::BasisFunctionType &basisfunction=initcond.basis.generator(initcond.e1);
    typename Basis1D::T y = initcond.R1pR2*x-initcond.R1;
    return    basisfunction(x,initcond.j1,initcond.k1,initcond.deriv1)
            * initcond.option.payoff_log(y)
            * exp(-2*initcond.eta*fabs(y))
            * initcond.sqrtR1pR2;
}

template <OptionType1D OType, typename Basis1D>
PayoffInitialCondition1D<OType,Basis1D>::PayoffInitialCondition1D(const Option1D<T,OType> &_option,
                                                                  const Basis1D &_basis,
                                                                  const T _eta,
                                                                  const T _R1, const T _R2)
    : option(_option), basis(_basis), eta(_eta), R1(_R1), R2(_R2), quadrature(*this),
      R1pR2(R1+R2), sqrtR1pR2(std::sqrt(R1+R2))
{

}

template <OptionType1D OType, typename Basis1D>
typename Basis1D::T
PayoffInitialCondition1D<OType,Basis1D>::operator()(XType _e1, int _j1, long _k1, int _deriv1) const
{
    j1 = _j1; k1 = _k1; e1 = _e1; deriv1 = _deriv1;
    return _integrate(*this);
}

template <OptionType1D OType, typename Basis1D>
typename Basis1D::T
PayoffInitialCondition1D<OType,Basis1D>::integrand(T x) const
{
    return _integrand(*this, x);
}

}   // namespace lawa
