namespace lawa {

template <OptionTypenD OType, ProcessType2D PType, typename Basis>
PayoffIntegral2D<OType,PType,Basis>::PayoffIntegral2D
(const Option2D<T,OType> &_option, ProcessParameters2D<T,PType> &_processparameters,
 const Basis &_basis, const T _left1, const T _right1, const T _left2, const T _right2)
    : option(_option), processparameters(_processparameters), basis(_basis),
      adapquad(*this),
      left1(_left1), right1(_right1), left2(_left2), right2(_right2), //quadrature(*this),
      RightmLeft1(right1-left1), SqrtRightmLeft1(std::sqrt(right1-left1)),
      RightmLeft2(right2-left2), SqrtRightmLeft2(std::sqrt(right2-left2))
{

}

template <OptionTypenD OType, ProcessType2D PType, typename Basis>
typename Basis::T
PayoffIntegral2D<OType,PType,Basis>::payoff(T x1, T x2) const
{
    return option.payoff_log(processparameters.u11*x1 + processparameters.u21*x2,
                             processparameters.u12*x1 + processparameters.u22*x2);
}

template <OptionTypenD OType, ProcessType2D PType, typename Basis>
typename Basis::T
PayoffIntegral2D<OType,PType,Basis>::integrand(T x1, T x2) const
{
    return 0.;
}

template <OptionTypenD OType, ProcessType2D PType, typename Basis>
typename Basis::T
PayoffIntegral2D<OType,PType,Basis>::operator()(Index2D &index2d)
{
    j1 = index2d.index1.j;
    j2 = index2d.index2.j;
    k1 = index2d.index1.k;
    k2 = index2d.index2.k;
    e1 = index2d.index1.xtype;
    e2 = index2d.index2.xtype;

    const typename Basis::FirstBasisType::BasisFunctionType  &first  = basis.first.generator(e1);
    const typename Basis::SecondBasisType::BasisFunctionType &second = basis.second.generator(e2);

    DenseVectorT singsupp1 = first.singularSupport(j1,k1);
    DenseVectorT singsupp2 = second.singularSupport(j2,k2);

    long double ret = 0.L;
    for (int i1=singsupp1.firstIndex(); i1<singsupp1.lastIndex(); ++i1) {
        T a1 = singsupp1(i1), b1 = singsupp1(i1+1);
        for (int i2=singsupp2.firstIndex(); i2<singsupp2.lastIndex(); ++i2) {
            T a2 = singsupp2(i2), b2 = singsupp2(i2+1);
            ret += (long double) adapquad.integrate(a1, b1, a2, b2);
        }
    }
    return (T)ret;
}

}   // namespace lawa
