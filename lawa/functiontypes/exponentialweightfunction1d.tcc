namespace lawa {

template <typename T>
T
ExponentialWeightFunction1D<T>::eta = 0.;

template <typename T>
T
ExponentialWeightFunction1D<T>::R1 = 0.;

template <typename T>
T
ExponentialWeightFunction1D<T>::R2 = 1.;

template <typename T>
flens::DenseVector<Array<T> >
ExponentialWeightFunction1D<T>::sing_pts(1);

template <typename T>
void
ExponentialWeightFunction1D<T>::setParameters(T _eta, T _R1, T _R2)
{
    //std::cout << "sing_pts = " << sing_pts << std::endl;
    eta=_eta;
    R1=_R1;
    R2=_R2;
}

template <typename T>
T
ExponentialWeightFunction1D<T>::weight(T x)
{
    return std::exp(-2*eta*fabs((R1+R2)*x-R1));
}

template <typename T>
T
ExponentialWeightFunction1D<T>::dweight(T x)
{
    T y = (R1+R2)*x-R1;
    if (y>0) {
        return -2.*eta*std::exp(-2*eta*fabs(y));
    }
    else {
        return  2.*eta*std::exp(-2*eta*fabs(y));
    }

}


}   // namespace lawa
