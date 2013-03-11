namespace lawa {

template<typename T, size_t PDim>
RB_System<T,PDim>::RB_System(ThetaStructure<T,PDim>& _thetas_a,
							 ThetaStructure<T,PDim>& _thetas_f)
 : thetas_a(_thetas_a), thetas_f(_thetas_f)
{}

template<typename T, size_t PDim>
void
RB_System<T,PDim>::
set_current_param(const std::array<T, PDim>& _param)
{
	current_param = _param;
}

template<typename T, size_t PDim>
std::array<T, PDim>&
RB_System<T,PDim>::
get_current_param()
{
	return current_param;
}
} // namespace lawa
