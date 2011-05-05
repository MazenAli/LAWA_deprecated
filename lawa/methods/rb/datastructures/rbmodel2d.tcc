namespace  lawa {

template <typename T, typename TruthModel>
RBModel2D<T, TruthModel>::RBModel2D()
{}


template <typename T, typename TruthModel>
void
RBModel2D<T, TruthModel>::set_current_param(const std::vector<T>& _param)
{
	current_param = _param;
}

template <typename T, typename TruthModel>
std::vector<T>&
RBModel2D<T, TruthModel>::get_current_param()
{
	return current_param;
}

template <typename T, typename TruthModel>
unsigned int
RBModel2D<T, TruthModel>::Q_a()
{
	return theta_a.size();	
}
        
template <typename T, typename TruthModel>
unsigned int
RBModel2D<T, TruthModel>::Q_f()
{
	return theta_f.size();	
}

template <typename T, typename TruthModel>
unsigned int
RBModel2D<T, TruthModel>::n_bf()
{
	return rb_basis_functions.size();	
}

template <typename T, typename TruthModel>
void
RBModel2D<T, TruthModel>::attach_inner_product_op(Operator2D<T>& _inner_product_op)
{
	inner_product_op = &_inner_product_op;
}

template <typename T, typename TruthModel>
void
RBModel2D<T, TruthModel>::set_truthmodel(TruthModel& _truthmodel)
{
	truth = &_truthmodel;
    truth->set_rb_model(*this);
}

template <typename T, typename TruthModel>
void
RBModel2D<T, TruthModel>::add_to_basis(const CoeffVector& sol)
{
	CoeffVector new_bf = sol;
    
    typename std::vector<CoeffVector>::iterator it;
    for (it = rb_basis_functions.begin(); it != rb_basis_functions.end(); ++it) {
        new_bf = new_bf - (*it) * inner_product((*it), sol); 
    }
    new_bf.scale(1./std::sqrt(inner_product(new_bf, new_bf)));
    
    rb_basis_functions.push_back(new_bf);
}

/*template <typename T, typename TruthSolver>
void
RBModel2D<T, TruthSolver>::update_RB_A_matrices()
{
	
}

template <typename T, typename TruthSolver>
void
RBModel2D<T, TruthSolver>::update_RB_F_vectors()
{	
    typename std::vector<CoeffVector>::const_iterator it;
    for (unsigned int q_f = 0; q_f <= Q_f(); q_f++) {
    	if (n_bf() == 1) {
    		RB_F_vectors[q_f] = new DenseVectorT(1);
        }
        else {
            DenseVectorT tmp(RB_F_vectors[q_f]);
            RB_F_vectors[q_f].engine().resize(n_bf());
            RB_F_vectors[q_f](tmp.range()) = tmp;
            RB_F_vectors[q_f](n_bf()) = 0.;            
        }

        for (it = rb_basis_functions[n_bf()-1].begin(); it != rb_basis_functions[n_bf()-1].end(); it++) {
            RB_F_vectors[q_f](n_bf()) += (*it).second * (*F_operators[q_f])((*it).first);
        }
    }
}*/

template <typename T, typename TruthModel>
T
RBModel2D<T, TruthModel>::inner_product(const CoeffVector& v1, const CoeffVector& v2)
{
	T val = 0;
    typename CoeffVector::const_iterator it1, it2;
    for (it1 = v1.begin(); it1 != v1.end() ; ++it1) {
        for (it2 = v2.begin(); it2 != v2.end(); ++it2) {
            val += (*it1).second * (*inner_product_op)((*it1).first, (*it2).first) * (*it2).second;
        }
    }
    return val;
}

} // namespace lawa

