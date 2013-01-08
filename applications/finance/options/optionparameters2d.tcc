namespace lawa {

template <typename T>
OptionParameters2D<T,BasketPut>::OptionParameters2D(void)
: strike(0.), maturity(0.), weight1(0.), weight2(0.),
  earlyExercise(false)
{

}

template <typename T>
OptionParameters2D<T,BasketPut>::OptionParameters2D
(T _strike, T _maturity, T _weight1, T _weight2, bool _earlyExercise)
: strike(_strike), maturity(_maturity), weight1(_weight1), weight2(_weight2),
  earlyExercise(_earlyExercise)
{

}

}   // namespace lawa
