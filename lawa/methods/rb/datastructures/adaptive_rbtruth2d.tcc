namespace  lawa {

template <typename T, typename Basis, typename TruthSolver>
AdaptiveRBTruth2D<T, Basis, TruthSolver>::AdaptiveRBTruth2D(Basis& _basis)
    : basis(_basis), lhs_op(this), rhs_op(this), repr_lhs_op(this), repr_rhs_A_op(this), repr_rhs_F_op(this)
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

template <typename T, typename Basis, typename TruthSolver>
RBModel2D<T, AdaptiveRBTruth2D<T, Basis, TruthSolver> >&
AdaptiveRBTruth2D<T, Basis, TruthSolver>::get_rb_model()
{
    return *rb;
}

template <typename T, typename Basis, typename TruthSolver>
void
AdaptiveRBTruth2D<T, Basis, TruthSolver>::calculate_representors()
{
    F_representors.resize(rb->Q_f());
    for (unsigned int i = 0; i < rb->Q_f(); ++i) {
        repr_rhs_F_op.set_current_op(*F_operators[i]);
        std::cout<< std::endl << " ==== Solving for Riesz Representor of F_" << i+1 << " =====" << std::endl<< std::endl;
        CoeffVector c = solver->repr_solve_F();
        F_representors[i] = c;
    }
    
    A_representors.resize(rb->n_bf());
    for (unsigned int n = 0; n < rb->n_bf(); ++n){
        repr_rhs_A_op.set_current_bf(rb->rb_basis_functions[n]);
        A_representors[n].resize(rb->Q_a());
        for (unsigned int i = 0; i < rb->Q_a(); ++i) {
            repr_rhs_A_op.set_current_op(*A_operators[i]);
            std::cout << std::endl<< " ==== Solving for Riesz Representor of A_" << i+1 
                      << "(" << n << ")" << " =====" << std::endl<< std::endl;
            A_representors[n][i] = solver->repr_solve_A();
        }
    }
}

template <typename T, typename Basis, typename TruthSolver>
void
AdaptiveRBTruth2D<T, Basis, TruthSolver>::calculate_representor_norms()
{
    // F_F representor norms 
    int Qf = rb->Q_f();
    rb->F_F_representor_norms.engine().resize(Qf, Qf);
    for(int qf1 = 1; qf1 <= Qf; ++qf1) {
        for (int qf2 = qf1; qf2 <= Qf; ++qf2) {
            rb->F_F_representor_norms(qf1,qf2) = rb->inner_product(F_representors[qf1-1], F_representors[qf2-1]);
            std::cout << " (F,F) = " << rb->F_F_representor_norms(qf1,qf2) << std::endl;
            if(qf1 != qf2) {
                rb->F_F_representor_norms(qf2,qf1) = rb->F_F_representor_norms(qf1,qf2);
            }
        }
    }
    
    // A_A representor norms 
    int Qa = rb->Q_a();
    int N = rb->n_bf();
    rb->A_A_representor_norms.resize(N);
    for(int n1 = 0; n1 < N; ++n1) {
        for(int n2 = 0; n2 < N; ++n2) {
            flens::GeMatrix<flens::FullStorage<T, cxxblas::ColMajor> > A_n1_n2(Qa, Qa);
            rb->A_A_representor_norms[n1].push_back(A_n1_n2);
            for(int qa1 = 1; qa1 <= Qa; ++qa1) {
                for(int qa2 = qa1; qa2 <= Qa; ++qa2) {
                    rb->A_A_representor_norms[n1][n2](qa1,qa2) 
                        = rb->inner_product(A_representors[n1][qa1-1], A_representors[n2][qa2-1]);
                    std::cout << " (A,A) = " << rb->A_A_representor_norms[n1][n2](qa1,qa2) << std::endl;
                    
                    if(qa1 != qa2) {
                        rb->A_A_representor_norms[n1][n2](qa2, qa1) = rb->A_A_representor_norms[n1][n2](qa1,qa2);
                    }
                }
            }
        }
    }
    
    // A_F representor norms 
    for(int n = 0; n < N; ++n) {
        flens::GeMatrix<flens::FullStorage<T, cxxblas::ColMajor> > A_F(Qa, Qf);
        rb->A_F_representor_norms.push_back(A_F);
        for(int qa = 1; qa <= Qa; ++qa) {
            for(int qf = 1; qf <= Qf; ++qf) {
                rb->A_F_representor_norms[n](qa, qf) = rb->inner_product(A_representors[n][qa-1], F_representors[qf-1]);
                std::cout << " (A,F) = " << rb->A_F_representor_norms[n](qa, qf) << std::endl;
            }
        }
    }
    
}

// ================================================================================================================ //
// ================================================================================================================ //

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

/*  Operator LHS_Representor */

template <typename T, typename Basis, typename TruthSolver>
T
AdaptiveRBTruth2D<T, Basis, TruthSolver>::Operator_LHS_Representor::
operator()(const Index2D &row_index, const Index2D &col_index)
{
    return (*thisTruth->rb->inner_product_op)(row_index, col_index);
}

/* Operator RHS_BilFormRepresentor */

template <typename T, typename Basis, typename TruthSolver>
T
AdaptiveRBTruth2D<T, Basis, TruthSolver>::Operator_RHS_BilFormRepresentor::operator()(const Index2D &lambda)
{
    T val = 0;
    typename Coefficients<Lexicographical,T,Index2D>::const_iterator it;
    for (it = current_bf->begin(); it != current_bf->end(); ++it) {
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
    typedef typename Coefficients<Lexicographical,T,Index2D>::value_type val_type;
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

template <typename T, typename Basis, typename TruthSolver>
void
AdaptiveRBTruth2D<T, Basis, TruthSolver>::Operator_RHS_BilFormRepresentor::
set_current_op(Operator2D<T>& op)
{
    current_op = &op;
}

template <typename T, typename Basis, typename TruthSolver>
void
AdaptiveRBTruth2D<T, Basis, TruthSolver>::Operator_RHS_BilFormRepresentor::
set_current_bf(Coefficients<Lexicographical,T,Index2D>& bf)
{
    current_bf = &bf;
}

/* Operator RHS_FunctionalRepresentor */

template <typename T, typename Basis, typename TruthSolver>
T
AdaptiveRBTruth2D<T, Basis, TruthSolver>::Operator_RHS_FunctionalRepresentor::operator()(const Index2D &lambda)
{
    return (*current_op)(lambda);
}

template <typename T, typename Basis, typename TruthSolver>
Coefficients<Lexicographical,T,Index2D>
AdaptiveRBTruth2D<T, Basis, TruthSolver>::Operator_RHS_FunctionalRepresentor::operator()(const IndexSet<Index2D> &Lambda)
{
    return (*current_op)(Lambda);
}

template <typename T, typename Basis, typename TruthSolver>
Coefficients<Lexicographical,T,Index2D>
AdaptiveRBTruth2D<T, Basis, TruthSolver>::Operator_RHS_FunctionalRepresentor::operator()(T tol)
{
    return (*current_op)(tol);
}

template <typename T, typename Basis, typename TruthSolver>
void
AdaptiveRBTruth2D<T, Basis, TruthSolver>::Operator_RHS_FunctionalRepresentor::
set_current_op(AdaptiveRhs<T, Index2D>& op)
{
    current_op = &op;
}

} // namespace lawa
