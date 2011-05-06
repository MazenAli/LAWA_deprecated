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
    
    update_RB_A_matrices();
    update_RB_F_vectors();
}

template <typename T, typename TruthSolver>
void
RBModel2D<T, TruthSolver>::update_RB_A_matrices()
{
	typename CoeffVector::const_iterator it, it1, it2;
    if ((n_bf() == 1) && (RB_A_matrices.size() < Q_a())) {
        //RB_A_matrices.resize(Q_a());
        for (unsigned int q_a = 0; q_a < Q_a(); ++q_a) {
            FullColMatrixT A(1,1);
            RB_A_matrices.push_back(A);
        }
    }
    for (unsigned int q_a = 0; q_a < Q_a(); ++q_a) {
        if (n_bf() == 1) {
            RB_A_matrices[q_a].engine().resize(1, 1);
        }
        else {
            FullColMatrixT tmp(RB_A_matrices[q_a]);
            RB_A_matrices[q_a].engine().resize((int)n_bf(), (int)n_bf());
            RB_A_matrices[q_a](tmp.rows(), tmp.cols()) = tmp;
            for(unsigned int i = 1; i <= n_bf(); ++i) {
                RB_A_matrices[q_a](i, n_bf()) = 0.;
                RB_A_matrices[q_a](n_bf(), i) = 0.;
            }
        }
        
        for (unsigned int i = 1; i < n_bf(); ++i) {
            for (it1 = rb_basis_functions[n_bf()-1].begin(); it1 != rb_basis_functions[n_bf()-1].end(); ++it1) {
                for (it2 = rb_basis_functions[i-1].begin(); it2 != rb_basis_functions[i-1].end(); ++it2) {
                    RB_A_matrices[q_a](n_bf(), i) += (*it1).second * (*it2).second 
                                                   * (*truth->A_operators[q_a])((*it1).first, (*it2).first);
                    RB_A_matrices[q_a](i, n_bf()) += (*it1).second * (*it2).second 
                                                   * (*truth->A_operators[q_a])((*it2).first, (*it1).first);
                }
            }        
        }
        for (it = rb_basis_functions[n_bf()-1].begin(); it != rb_basis_functions[n_bf()-1].end(); ++it) {
            RB_A_matrices[q_a](n_bf(), n_bf()) += (*it).second * (*it).second 
                                                * (*truth->A_operators[q_a])((*it).first, (*it).first);            
        }

    }
}

template <typename T, typename TruthSolver>
void
RBModel2D<T, TruthSolver>::update_RB_F_vectors()
{	
    typename CoeffVector::const_iterator it;
    if ((n_bf() == 1) && (RB_F_vectors.size() < Q_f())) {
        RB_F_vectors.resize(Q_f());
    }
    for (unsigned int q_f = 0; q_f < Q_f(); ++q_f) {
    	if (n_bf() == 1) {
    		RB_F_vectors[q_f].engine().resize(1);
        }
        else {
            DenseVectorT tmp(RB_F_vectors[q_f]);
            RB_F_vectors[q_f].engine().resize((int)n_bf());
            RB_F_vectors[q_f](tmp.range()) = tmp;
            RB_F_vectors[q_f](n_bf()) = 0.;            
        }

        for (it = rb_basis_functions[n_bf()-1].begin(); it != rb_basis_functions[n_bf()-1].end(); ++it) {
            RB_F_vectors[q_f](n_bf()) += (*it).second * (*truth->F_operators[q_f])((*it).first);
        }
    }
}

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

