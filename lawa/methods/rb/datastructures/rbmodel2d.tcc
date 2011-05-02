namespace  lawa {

template <typename T, typename TruthSolver>
RBModel2D<T, TruthSolver>::RBModel2D()
{}


template <typename T, typename TruthSolver>
void
RBModel2D<T, TruthSolver>::set_current_param(const std::vector<T>& _param)
{
	current_param = _param;
}

template <typename T, typename TruthSolver>
std::vector<T>&
RBModel2D<T, TruthSolver>::get_current_param()
{
	return current_param;
}

} // namespace lawa

