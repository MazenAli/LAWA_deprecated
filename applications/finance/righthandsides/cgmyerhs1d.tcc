namespace lawa {


template <typename T, typename Basis1D>
OptionRHS1D<T, Put, CGMYe, Basis1D>::OptionRHS1D(const OptionParameters1D<T,Put> &_optionparameters,
                                                const ProcessParameters1D<T,CGMYe> &_processparameters,
                                                const Basis1D &_basis, T _R1, T _R2,
                                                bool _excessToPayoff)
: optionparameters(_optionparameters), processparameters(_processparameters),
  basis(_basis), kernel(_processparameters),
  R1(_R1), R2(_R2),
  OneDivSqrtR2pR1(1./std::sqrt(R2+R1)), OneDivR2pR1(1./(R2+R1)), R1DivR1pR2(R1/(R1+R2)),
  excessToPayoff(_excessToPayoff)
{
    if (R1==0 && R2==1) {
        OneDivSqrtR2pR1=1.;
        OneDivR2pR1    =1.;
        R1DivR1pR2     =0.;
    }
    bool is_realline_basis = flens::IsSame<Basis<T,Primal,R,CDF>, Basis1D>::value;
    assert(is_realline_basis || (R1>0 && R2>0.));
}

template <typename T, typename Basis1D>
T
OptionRHS1D<T, Put, CGMYe, Basis1D>::operator()(XType xtype, int j, int k) const
{
    GeMatrix<FullStorage<T,ColMajor> > varphi_mu_deltas;
    T K = optionparameters.strike;
    T sigma = processparameters.sigma;

    if (!excessToPayoff) return 0.;

    if (basis.d==2) {
        T ret = 0.;
        T step = pow2i<T>(-j-5);
        ret += K*(kernel.ExpXmOne_k1_pos + kernel.ExpXmOne_k1_neg)
                  * OneDivSqrtR2pR1 * basis.generator(xtype)(R1DivR1pR2,j,k,0);
        ret -= K*kernel.ExpXmOne_k2_neg
                * OneDivR2pR1*OneDivSqrtR2pR1 * basis.generator(xtype)(R1DivR1pR2+step,j,k,1);
        ret -= K*kernel.ExpXmOne_k2_pos
                * OneDivR2pR1*OneDivSqrtR2pR1 * basis.generator(xtype)(R1DivR1pR2-step,j,k,1);

        GeMatrix<FullStorage<T,ColMajor> > varphi_row_deltas;
        varphi_row_deltas = computeDeltas<T,Basis1D>(basis,j,k,xtype);

        if (R1!=0 && R2!=1) {
            varphi_row_deltas(_,1)  *=(R1+R2);  varphi_row_deltas(_,1)-=R1;
            varphi_row_deltas(_,2)  *= OneDivR2pR1*OneDivSqrtR2pR1;
        }

        for (int i=varphi_row_deltas.rows().firstIndex(); i<=varphi_row_deltas.rows().lastIndex(); ++i) {
            T x = varphi_row_deltas(i,1);
            T c = varphi_row_deltas(i,2);
            if (x != 0) {
                ret += c*( K*exp(x)*kernel.ThirdTailFirstExpMomIntegral(-x)
                       - K*kernel.ThirdTailIntegral(-x));
            }
        }

        ret += 0.5*sigma*sigma* K * OneDivSqrtR2pR1 * basis.generator(xtype)(R1DivR1pR2,j,k,0);

        return ret;
    }
    else {
        assert(0);
        return 0;
    }
}

template <typename T, typename Basis1D>
T
OptionRHS1D<T, Put, CGMYe, Basis1D>::operator()(T time, XType xtype, int j, int k) const
{
    return this->operator()(xtype, j, k);
}

template <typename T, typename Basis1D>
T
OptionRHS1D<T, Put, CGMYe, Basis1D>::operator()(const Index1D &lambda) const
{
    return this->operator()(lambda.xtype, lambda.j, lambda.k);
}

template <typename T, typename Basis1D>
T
OptionRHS1D<T, Put, CGMYe, Basis1D>::operator()(T time, const Index1D &lambda) const
{
    return this->operator()(lambda.xtype, lambda.j, lambda.k);
}



}   // namespace lawa

