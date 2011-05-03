namespace  lawa {

template <typename T, typename Basis, typename TruthSolver>
AdaptiveRBModel2D<T, Basis, TruthSolver>::AdaptiveRBModel2D(Basis& _basis)
	: RBModel2D<T, TruthSolver>(), basis(_basis), lhs_op(this), rhs_op(this)
{}

template <typename T, typename Basis, typename TruthSolver>
void
AdaptiveRBModel2D<T, Basis, TruthSolver>::attach_A_q(theta_fctptr theta_a_q, Operator2D<T>& A_q)
{
	this->theta_a.push_back(theta_a_q);
	A_operators.push_back(&A_q);
}

template <typename T, typename Basis, typename TruthSolver>
void
AdaptiveRBModel2D<T, Basis, TruthSolver>::attach_F_q(theta_fctptr theta_f_q, AdaptiveRhs<T, Index2D>& F_q)
{
	this->theta_f.push_back(theta_f_q);
	F_operators.push_back(&F_q);
}

template <typename T, typename Basis, typename TruthSolver>
void 
AdaptiveRBModel2D<T, Basis, TruthSolver>::set_truthsolver(TruthSolver& _truthsolver)
{
	this->truthsolver = &_truthsolver;
    this->truthsolver->set_model(*this);
}


/*  Operator LHS */

template <typename T, typename Basis, typename TruthSolver>
T
AdaptiveRBModel2D<T, Basis, TruthSolver>::Operator_LHS::operator()(const Index2D &row_index, const Index2D &col_index)
{
	T val = 0;
	for (unsigned int i = 0; i < thisModel->A_operators.size(); ++i) {
        val += (*thisModel->theta_a[i])(thisModel->current_param) 
        	 * (*thisModel->A_operators[i])(row_index, col_index);
    }
    
    return val;
}

/*  Operator RHS */

template <typename T, typename Basis, typename TruthSolver>
T
AdaptiveRBModel2D<T, Basis, TruthSolver>::Operator_RHS::operator()(const Index2D &lambda)
{
	T val = 0;
	for (unsigned int i = 0; i < thisModel->F_operators.size(); ++i) {
        val += (*thisModel->theta_f[i])(thisModel->current_param) 
        	 * (*thisModel->F_operators[i])(lambda);
    }
    
    return val;
}

template <typename T, typename Basis, typename TruthSolver>
Coefficients<Lexicographical,T,Index2D>
AdaptiveRBModel2D<T, Basis, TruthSolver>::Operator_RHS::operator()(const IndexSet<Index2D> &Lambda)
{
	Coefficients<Lexicographical,T,Index2D> c;
	for (unsigned int i = 0; i < thisModel->F_operators.size(); ++i) {
        c = c + ((*thisModel->theta_f[i])(thisModel->current_param) * (*thisModel->F_operators[i])(Lambda));
    }
    
    return c;
}

template <typename T, typename Basis, typename TruthSolver>
Coefficients<Lexicographical,T,Index2D>
AdaptiveRBModel2D<T, Basis, TruthSolver>::Operator_RHS::operator()(T tol)
{
	Coefficients<Lexicographical,T,Index2D> c;
	for (unsigned int i = 0; i < thisModel->F_operators.size(); ++i) {
        c += (*thisModel->theta_f[i])(thisModel->current_param) 
        	 * (*thisModel->F_operators[i])(tol);
    }
    
    return c;
}

} // namespace lawa
