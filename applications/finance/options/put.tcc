namespace lawa {

template <typename T>
Option1D<T,Put>::Option1D(const OptionParameters1D<T,Put> &_optionparameters)
    : optionparameters(_optionparameters)

{
    assert(optionparameters.strike>=0.);
    assert(optionparameters.maturity>=0.);
    singularPoints.engine().resize(1);
    singularPoints(1) = 0.;
}

template <typename T>
T
Option1D<T,Put>::payoff(T S) const
{
    return std::max( optionparameters.strike - S, 0.0);
}

template <typename T>
T
Option1D<T,Put>::payoff_log(T x) const
{
    return optionparameters.strike*std::max(1.-std::exp(x), 0.);
}


template <typename T>
T
Option1D<T,Put>::value(const ProcessParameters1D<T,BlackScholes> &processparameters, T S, T t)
{
    assert(optionparameters.earlyExercise==false);
    static boost::math::normal norm;
    T r        = processparameters.r;
    T sigma    = processparameters.sigma;
    T strike   = optionparameters.strike;
    T maturity = optionparameters.maturity;
    T d_1 = ( log(S/strike) + ( r + sigma*sigma/2.0 ) * (maturity-t) )
                 / ( sigma * sqrt(maturity-t) );
    T d_2 = d_1 - sigma * sqrt(maturity-t);
    return strike*exp(-r*(maturity-t))*cdf(norm, -d_2) - S * cdf(norm, -d_1);
}

template <typename T>
T
Option1D<T,Put>::value(const ProcessParameters1D<T,CGMY> &processparameters, T S, T t)
{
    assert(optionparameters.earlyExercise==false);
    assert(processparameters.k_Y!=0);
    assert(processparameters.k_Y!=1);
    CharacteristicFunction1D<T,CGMY> charfunc(processparameters);

    FourierPricer1D<T, CGMY> fp(charfunc, S, optionparameters.maturity,
                                std::max(optionparameters.strike-10.,0.),
                                optionparameters.strike+10.);
    fp.solve(40000,17);
    return fp(optionparameters.strike);
}

}   // namespace lawa