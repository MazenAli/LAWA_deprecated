namespace lawa {

template <typename T>
OptionParameters1D<T,Put>::OptionParameters1D(T _strike, T _maturity, T _earlyExercise)
: strike(_strike), maturity(_maturity), earlyExercise(_earlyExercise)
{

}

}   // namespace lawa
