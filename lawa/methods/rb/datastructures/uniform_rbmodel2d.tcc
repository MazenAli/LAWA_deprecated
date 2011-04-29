#include <lawa/methods/uniform/algorithms/assembler2d.h>

namespace  lawa {

template <typename T, typename Prec>
UniformRBModel2D<T, Prec>::UniformRBModel2D()
	: RBModel2D<T>(), lhs_op(this), rhs_op(this), noprec(), prec(noprec)
{}

template <typename T, typename Prec>
UniformRBModel2D<T, Prec>::UniformRBModel2D(Prec& _prec)
	: RBModel2D<T>(), lhs_op(this), rhs_op(this), noprec(), prec(_prec)
{}


template <typename T, typename Prec>
void
UniformRBModel2D<T, Prec>::attach_A_q(theta_fctptr theta_a_q, Operator2D<T>& A_q)
{
	this->theta_a.push_back(theta_a_q);
	A_operators.push_back(&A_q);
}


template <typename T, typename Prec>
void
UniformRBModel2D<T, Prec>::attach_F_q(theta_fctptr theta_f_q, Rhs2D<T>& F_q)
{
	this->theta_f.push_back(theta_f_q);
	F_operators.push_back(&F_q);
}


template <typename T, typename Prec>
T
UniformRBModel2D<T, Prec>::Operator_LHS::operator()(XType row_xtype_x, int j1_x, int k1_x,
                                                    XType row_xtype_y, int j1_y, int k1_y,
                                       		        XType col_xtype_x, int j2_x, int k2_x,
                                         	        XType col_xtype_y, int j2_y, int k2_y) const
{
	T val = 0;
	for (unsigned int i = 0; i < thisModel->A_operators.size(); ++i) {
        val += (*thisModel->theta_a[i])(thisModel->current_param) 
        	 * (*thisModel->A_operators[i])(row_xtype_x, j1_x, k1_x, row_xtype_y, j1_y, k1_y,
                                 col_xtype_x, j2_x, k2_x, col_xtype_y, j2_y, k2_y);
    }
    
    return val;
}

template <typename T, typename Prec>
T
UniformRBModel2D<T, Prec>::Operator_RHS::operator()(XType xtype_x, int j_x, int k_x,
                                         	        XType xtype_y, int j_y, int k_y) const
{
	T val = 0;
	for (unsigned int i = 0; i < thisModel->F_operators.size(); ++i) {
        val += (*thisModel->theta_f[i])(thisModel->current_param) 
        	 * (*thisModel->F_operators[i])(xtype_x, j_x, k_x, xtype_y, j_y, k_y);
    }
    
    return val;
}

} // namespace lawa
