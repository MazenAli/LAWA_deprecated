namespace  lawa {

template <typename T, typename Basis, typename TruthSolver>
AdaptiveRBTruth2D<T, Basis, TruthSolver>::AdaptiveRBTruth2D(Basis& _basis)
    : basis(_basis), lhs_op(this), rhs_op(this)
{}

template <typename T, typename Basis, typename TruthSolver>
void
AdaptiveRBTruth2D<T, Basis, TruthSolver>::attach_A_q(theta_fctptr theta_a_q, Operator2D<T>& A_q)
{
    rb->theta_a.push_back(theta_a_q);
    A_operators.push_back(&A_q);
}

template <typename T, typename Basis, typename TruthSolver>
void
AdaptiveRBTruth2D<T, Basis, TruthSolver>::attach_F_q(theta_fctptr theta_f_q, AdaptiveRhs<T, Index2D>& F_q)
{
    rb->theta_f.push_back(theta_f_q);
    F_operators.push_back(&F_q);
}

template <typename T, typename Basis, typename TruthSolver>
void 
AdaptiveRBTruth2D<T, Basis, TruthSolver>::set_truthsolver(TruthSolver& _truthsolver)
{
    solver = &_truthsolver;
    solver->set_model(*this);
}

template <typename T, typename Basis, typename TruthSolver>
void
AdaptiveRBTruth2D<T, Basis, TruthSolver>::set_rb_model(RBModel2D<T, AdaptiveRBTruth2D<T, Basis, TruthSolver> >& _rb)
{
    rb = &_rb;
}


/*  Operator LHS */

template <typename T, typename Basis, typename TruthSolver>
T
AdaptiveRBTruth2D<T, Basis, TruthSolver>::Operator_LHS::operator()(const Index2D &row_index, const Index2D &col_index)
{
    T val = 0;
    for (unsigned int i = 0; i < thisTruth->A_operators.size(); ++i) {
        val += (*thisTruth->rb->theta_a[i])(thisTruth->rb->get_current_param()) 
             * (*thisTruth->A_operators[i])(row_index, col_index);
    }
    
    return val;
}

/*  Operator RHS */

template <typename T, typename Basis, typename TruthSolver>
T
AdaptiveRBTruth2D<T, Basis, TruthSolver>::Operator_RHS::operator()(const Index2D &lambda)
{
    T val = 0;
    for (unsigned int i = 0; i < thisTruth->F_operators.size(); ++i) {
        val += (*thisTruth->rb->theta_f[i])(thisTruth->rb->get_current_param()) 
             * (*thisTruth->F_operators[i])(lambda);
    }
    
    return val;
}

template <typename T, typename Basis, typename TruthSolver>
Coefficients<Lexicographical,T,Index2D>
AdaptiveRBTruth2D<T, Basis, TruthSolver>::Operator_RHS::operator()(const IndexSet<Index2D> &Lambda)
{
    Coefficients<Lexicographical,T,Index2D> c;
    for (unsigned int i = 0; i < thisTruth->F_operators.size(); ++i) {
        c = c + (*thisTruth->F_operators[i])(Lambda) 
              * (*thisTruth->rb->theta_f[i])(thisTruth->rb->get_current_param());
    }
    
    return c;
}

template <typename T, typename Basis, typename TruthSolver>
Coefficients<Lexicographical,T,Index2D>
AdaptiveRBTruth2D<T, Basis, TruthSolver>::Operator_RHS::operator()(T tol)
{
    Coefficients<Lexicographical,T,Index2D> c;
    for (unsigned int i = 0; i < thisTruth->F_operators.size(); ++i) {
        c += (*thisTruth->F_operators[i])(tol) 
           * (*thisTruth->rb->theta_f[i])(thisTruth->rb->get_current_param());
    }
    
    return c;
}

/* Operator RHS_BilFormRepresentor */

template <typename T, typename Basis, typename TruthSolver>
T
AdaptiveRBTruth2D<T, Basis, TruthSolver>::Operator_RHS_BilFormRepresentor::operator()(const Index2D &lambda)
{
    T val = 0;
    typename Coefficients<Lexicographical,T,Index2D>::const_iterator it;
    for (it = current_bf.begin(); it != current_bf.end(); ++it) {
        val += (*it).second * (*current_op)((*it).first, lambda);
    }
    
    return - val;
}

template <typename T, typename Basis, typename TruthSolver>
Coefficients<Lexicographical,T,Index2D>
AdaptiveRBTruth2D<T, Basis, TruthSolver>::Operator_RHS_BilFormRepresentor::operator()(const IndexSet<Index2D> &Lambda)
{
    Coefficients<Lexicographical, T, Index2D> coeffs;
    
    typename IndexSet<Index2D>::const_iterator it_Lambda;
    typedef typename Coefficients<Lexicographical,T,Index>::value_type val_type;
    for (it_Lambda = Lambda.begin(); it_Lambda != Lambda.end(); ++it_Lambda) {
        coeffs.insert(val_type((*it_Lambda), this->operator()(*it_Lambda)));
    }
    
    return coeffs;
}

template <typename T, typename Basis, typename TruthSolver>
Coefficients<Lexicographical,T,Index2D>
AdaptiveRBTruth2D<T, Basis, TruthSolver>::Operator_RHS_BilFormRepresentor::operator()(T tol)
{
    Coefficients<Lexicographical, T, Index2D> coeffs;
    
    std::cerr << "Operator()(T tol) : Not implemented yet" << std::endl;
    exit(1);
    
    return coeffs;
}

} // namespace lawa
