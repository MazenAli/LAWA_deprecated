namespace lawa {

template <typename T>
Option2D<T,SumOfPuts>::Option2D(void)
    : optionparameters()

{

}

template <typename T>
Option2D<T,SumOfPuts>::Option2D(const OptionParameters2D<T,SumOfPuts> &_optionparameters)
    : optionparameters(_optionparameters)

{
    assert(optionparameters.strike1>0.);
    assert(optionparameters.strike2>0.);
    assert(optionparameters.maturity>=0.);
}

template <typename T>
T
Option2D<T,SumOfPuts>::payoff(T S1, T S2) const
{
    return    optionparameters.w1*std::max(optionparameters.strike1 - S1)
            + optionparameters.w2*std::max(optionparameters.strike2 - S2);
}

template <typename T>
T
Option2D<T,SumOfPuts>::payoff_log(T x1, T x2) const
{
    return   optionparameters.weight1 * optionparameters.strike1*std::max(1.-std::exp(x1),(T)0.)
           + optionparameters.weight2 * optionparameters.strike2*std::max(1.-std::exp(x2),(T)0.);
}

}   //namespace lawa
