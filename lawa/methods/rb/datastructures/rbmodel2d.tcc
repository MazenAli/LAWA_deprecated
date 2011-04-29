namespace  lawa {

template <typename T>
RBModel2D<T>::RBModel2D()
{}


template <typename T>
void
RBModel2D<T>::set_current_param(const std::vector<T>& _param)
{
	current_param = _param;
}

template <typename T>
std::vector<T>&
RBModel2D<T>::get_current_param()
{
	return current_param;
}

} // namespace lawa

