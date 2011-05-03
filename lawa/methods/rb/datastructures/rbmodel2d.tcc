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

template <typename T, typename TruthSolver>
void
RBModel2D<T, TruthSolver>::attach_inner_product_op(Operator2D<T>& _inner_product_op)
{
	inner_product_op = &_inner_product_op;
}

template <typename T, typename TruthSolver>
void
RBModel2D<T, TruthSolver>::add_to_basis(CoeffVector& sol)
{
	CoeffVector new_bf = sol;
    
    std::vector<CoeffVector>::const_iterator it;
    for (it = rb_basis_functions.begin(); it != rb_basis_functions.end(); ++it) {
        sol = sol; 
    }
}

} // namespace lawa

