namespace  lawa {

template <typename T, typename Basis>
AdaptiveRBModel2D<T, Basis>::AdaptiveRBModel2D()
	: RBModel<T>()
{}

template <typename T, typename Basis>
void
AdaptiveRBModel2D<T, Basis>::attach_A_q(theta_fctptr theta_a_q, Operator2D<T>& A_q)
{
	this->theta_a.push_back(theta_a_q);
	A_operators.push_back(&A_q);
}

template <typename T, typename Basis>
void
AdaptiveRBModel2D<T, Basis>::attach_F_q(theta_fctptr theta_f_q, Rhs2D<T>& F_q)
{
	this->theta_f.push_back(theta_f_q);
	F_operators.push_back(&F_q);
}



} // namespace lawa