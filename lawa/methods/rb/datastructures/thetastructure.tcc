namespace lawa {

template<typename T,size_t PDim>
ThetaStructure<T,PDim>::ThetaStructure()
 : current_param(nullptr){}

template<typename T,size_t PDim>
ThetaStructure<T,PDim>::ThetaStructure(const std::vector<ThetaFct>& _thetas)
 : thetas(_thetas), current_param(nullptr)
{}


template<typename T,size_t PDim>
unsigned int
ThetaStructure<T,PDim>::size() const
{
	return thetas.size();
}

template<typename T,size_t PDim>
void
ThetaStructure<T,PDim>::set_current_param(const std::array<T, PDim>& _param)
{
	current_param = _param;
}


template<typename T,size_t PDim>
std::array<T, PDim>&
ThetaStructure<T,PDim>::get_current_param()
{
	return current_param;
}

template<typename T,size_t PDim>
T
ThetaStructure<T,PDim>::eval(int i, std::array<T,PDim>& mu)
{
	return (thetas[i])(mu);
}

template<typename T,size_t PDim>
T
ThetaStructure<T,PDim>::eval(int i)
{
	return eval(i,*current_param);
}

} // namespace lawa
