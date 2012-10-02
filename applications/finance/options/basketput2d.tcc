namespace lawa {

template <typename T>
Option2D<T,BasketPut>::Option2D(const OptionParameters2D<T,BasketPut> &_optionparameters)
    : optionparameters(_optionparameters)

{
    assert(optionparameters.strike>=0.);
    assert(optionparameters.maturity>=0.);
}

template <typename T>
T
Option2D<T,BasketPut>::payoff(T S1, T S2) const
{
    return std::max(optionparameters.strike - optionparameters.w1*S1 - optionparameters.w2*S2, 0.0);
}

template <typename T>
T
Option2D<T,BasketPut>::payoff_log(T x1, T x2) const
{
    return optionparameters.strike*std::max(1.-optionparameters.weight1*std::exp(x1)
                                              -optionparameters.weight2*std::exp(x2), 0.);
}

}   //namespace lawa
