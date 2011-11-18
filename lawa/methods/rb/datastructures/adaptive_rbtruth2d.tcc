#include <lawa/methods/rb/postprocessing/plotting.h>

namespace  lawa {

template <typename T, typename Basis, typename Prec, typename TruthSolver, typename Compression>
AdaptiveRBTruth2D<T, Basis, Prec, TruthSolver, Compression>::
AdaptiveRBTruth2D(Basis& _basis, Prec& _prec, bool _use_inner_product, bool _use_A_matrix)
    : basis(_basis), lhs_op(this), rhs_op(this), repr_lhs_op(this), repr_rhs_A_op(this), repr_rhs_F_op(this),
      use_inner_product_matrix(_use_inner_product), use_A_operator_matrices(_use_A_matrix), 
      assembled_inner_product_matrix(false), assembled_prec_vec(false), assembled_A_operator_matrices(false),
      prec(_prec), prec_data()
{}

template <typename T, typename Basis, typename Prec, typename TruthSolver, typename Compression>
void
AdaptiveRBTruth2D<T, Basis, Prec, TruthSolver, Compression>::attach_A_q(theta_fctptr theta_a_q, Operator2D<T>& A_q)
{
    rb->theta_a.push_back(theta_a_q);
    A_operators.push_back(&A_q);
}

template <typename T, typename Basis, typename Prec, typename TruthSolver, typename Compression>
void
AdaptiveRBTruth2D<T, Basis, Prec, TruthSolver, Compression>::attach_A_q(Operator2D<T>& A_q)
{
    A_operators.push_back(&A_q);
}

template <typename T, typename Basis, typename Prec, typename TruthSolver, typename Compression>
void
AdaptiveRBTruth2D<T, Basis, Prec, TruthSolver, Compression>::attach_F_q(theta_fctptr theta_f_q, AdaptiveRhs<T, Index2D>& F_q)
{
    rb->theta_f.push_back(theta_f_q);
    F_operators.push_back(&F_q);
}

template <typename T, typename Basis, typename Prec, typename TruthSolver, typename Compression>
void
AdaptiveRBTruth2D<T, Basis, Prec, TruthSolver, Compression>::attach_F_q(AdaptiveRhs<T, Index2D>& F_q)
{
    F_operators.push_back(&F_q);
}

template <typename T, typename Basis, typename Prec, typename TruthSolver, typename Compression>
void
AdaptiveRBTruth2D<T, Basis, Prec, TruthSolver, Compression>::attach_inner_product_op(Operator2D<T>& _inner_product_op)
{
    inner_product_op = &_inner_product_op;
}

template <typename T, typename Basis, typename Prec, typename TruthSolver, typename Compression>
void 
AdaptiveRBTruth2D<T, Basis, Prec, TruthSolver, Compression>::set_truthsolver(TruthSolver& _truthsolver)
{
    solver = &_truthsolver;
    solver->set_model(*this);
}

template <typename T, typename Basis, typename Prec, typename TruthSolver, typename Compression>
void
AdaptiveRBTruth2D<T, Basis, Prec, TruthSolver, Compression>::
set_rb_model(RBModel2D<T, AdaptiveRBTruth2D<T, Basis, Prec, TruthSolver, Compression> >& _rb)
{
    rb = &_rb;
}

template <typename T, typename Basis, typename Prec, typename TruthSolver, typename Compression>
RBModel2D<T, AdaptiveRBTruth2D<T, Basis, Prec, TruthSolver, Compression> >&
AdaptiveRBTruth2D<T, Basis, Prec, TruthSolver, Compression>::get_rb_model()
{
    return *rb;
}

template <typename T, typename Basis, typename Prec, typename TruthSolver, typename Compression>
Coefficients<Lexicographical,T,Index2D>
AdaptiveRBTruth2D<T, Basis, Prec, TruthSolver, Compression>::truth_solve()
{
    CoeffVector u = solver->truth_solve();
    
    undo_prec(u);
    
    return u;
}

template <typename T, typename Basis, typename Prec, typename TruthSolver, typename Compression>
void
AdaptiveRBTruth2D<T, Basis, Prec, TruthSolver, Compression>::undo_prec(CoeffVector& u)
{    
    // Undo preconditioning
    typename CoeffVector::iterator it;
    for(it = u.begin(); it != u.end(); ++it){
        (*it).second *= get_prec((*it).first);
    }     
}

template <typename T, typename Basis, typename Prec, typename TruthSolver, typename Compression>
T
AdaptiveRBTruth2D<T, Basis, Prec, TruthSolver, Compression>::get_prec(const Index2D& index)
{
    typedef typename Coefficients<Lexicographical,T,Index2D>::const_iterator const_coeff_it;
    T precval = 1.;
    
    if (!flens::IsSame<NoPreconditioner<T,Index2D>, Prec>::value) {
        // Right precondioning:
        const_coeff_it it_index   = prec_data.find(index);
        //  Entry has already been computed:
        if (it_index != prec_data.end()) {
            precval *= (*it_index).second;
        }
        //  Entry has not yet been computed:
        else {
            T tmp = prec(index);
            prec_data[index] = tmp;
            precval *= tmp;
        }
    }
    return precval;
}

template <typename T, typename Basis, typename Prec, typename TruthSolver, typename Compression>
void
AdaptiveRBTruth2D<T, Basis, Prec, TruthSolver, Compression>::add_new_basis_function(const CoeffVector& sol)
{
    typename std::vector<CoeffVector>::iterator it;
    CoeffVector new_bf = sol;

    CoeffVector bf;
    for (it = rb->rb_basis_functions.begin(); it != rb->rb_basis_functions.end(); ++it) {
        bf = (*it);
        new_bf = new_bf - bf * inner_product(bf, sol); 
    }
    
    std::cout << "Gram - Schmidt: norm^2 = " << inner_product(new_bf, new_bf) << std::endl;
    new_bf.scale(1./std::sqrt(inner_product(new_bf, new_bf)));
    rb->rb_basis_functions.push_back(new_bf);    
}

template <typename T, typename Basis, typename Prec, typename TruthSolver, typename Compression>
void
AdaptiveRBTruth2D<T, Basis, Prec, TruthSolver, Compression>::update_representors()
{
    int N = rb->n_bf();
    if(N == 1) {
        F_representors.resize(rb->Q_f());
        for (unsigned int i = 0; i < rb->Q_f(); ++i) {
            repr_rhs_F_op.set_current_op(*F_operators[i]);
            std::cout<< std::endl << "  --- Solving for Riesz Representor of F_" << i+1 << "---" << std::endl<< std::endl;
            CoeffVector c = solver->repr_solve_F();
            
            undo_prec(c);
            
            F_representors[i] = c;
        } 
    }
    
    A_representors.resize(N);
    repr_rhs_A_op.set_current_bf(rb->rb_basis_functions[N-1]);
    A_representors[N-1].resize(rb->Q_a());
    for (unsigned int i = 0; i < rb->Q_a(); ++i) {
        repr_rhs_A_op.set_current_op(*A_operators[i]);
        std::cout << std::endl<< "  --- Solving for Riesz Representor of A_" << i+1 
                  << "(" << N << ")" << " --- " << std::endl<< std::endl;
         
        CoeffVector c = solver->repr_solve_A();   
        
        undo_prec(c);
            
        A_representors[N-1][i] = c;
    }
 }

template <typename T, typename Basis, typename Prec, typename TruthSolver, typename Compression>
void
AdaptiveRBTruth2D<T, Basis, Prec, TruthSolver, Compression>::update_representor_norms()
{
    Timer timer;
    int N = rb->n_bf();
    int Qf = rb->Q_f();
    
    CoeffVector vec1, vec2;

    if(N == 1){
      std::cout << "  ... F x F .... " << std::endl;
       // F_F representor norms 
       rb->F_F_representor_norms.engine().resize(Qf, Qf);

       timer.start();
       for(int qf1 = 1; qf1 <= Qf; ++qf1) {
           for (int qf2 = qf1; qf2 <= Qf; ++qf2) {
               
               vec1 = F_representors[qf1-1];
               vec2 = F_representors[qf2-1];

               rb->F_F_representor_norms(qf1,qf2) = inner_product(vec1, vec2);

               if(qf1 != qf2) {
                   rb->F_F_representor_norms(qf2,qf1) = rb->F_F_representor_norms(qf1,qf2);
               }
           }
       }   
       timer.stop();
       std::cout << "  " << timer.elapsed() << " seconds" << std::endl << std::endl;
       std::cout << rb->F_F_representor_norms << std::endl;
    }
     
     
    std::cout << "  ...  A x A .... " << std::endl;    
     
     
    // A_A representor norms 
    int Qa = rb->Q_a();
    timer.start();
    for(int n1 = 0; n1 < N; ++n1) {
        flens::GeMatrix<flens::FullStorage<T, cxxblas::ColMajor> > A_n1_N(Qa, Qa);
        for(int qa1 = 1; qa1 <= Qa; ++qa1) {
          for(int qa2 = qa1; qa2 <= Qa; ++qa2) {  
          
            vec1 = A_representors[n1][qa1-1];
            vec2 = A_representors[N-1][qa2-1];
                              
            A_n1_N(qa1, qa2) = inner_product(vec1, vec2);
            
            if(qa1 != qa2){
              if(n1 == N-1){
                A_n1_N(qa2, qa1) = A_n1_N(qa1, qa2);
              }
              else{
              
                vec1 = A_representors[n1][qa2-1];
                vec2 = A_representors[N-1][qa1-1];
                  
                A_n1_N(qa2, qa1) = inner_product(vec1, vec2);                              
              }
            }
          }
        }
        if(n1 == N-1){
          std::vector<flens::GeMatrix<flens::FullStorage<T, cxxblas::ColMajor> > > newvec;
          newvec.push_back(A_n1_N);
          rb->A_A_representor_norms.push_back(newvec);
        }
        else{
          rb->A_A_representor_norms[n1].push_back(A_n1_N);          
        }
        std::cout << "n1 = " << n1 << ", n2 = " << N-1 << " " <<  rb->A_A_representor_norms[n1][N-1-n1] << std::endl;
    }
    timer.stop();
    std::cout << "  " << timer.elapsed() << " seconds" << std::endl << std::endl;

    // A_F representor norms 
    std::cout << "  ... A x F .... " << std::endl;
    timer.start();
    flens::GeMatrix<flens::FullStorage<T, cxxblas::ColMajor> > A_F(Qa, Qf);
    for(int qa = 1; qa <= Qa; ++qa) {
        for(int qf = 1; qf <= Qf; ++qf) {
           
            vec1 = A_representors[N-1][qa-1];
            vec2 = F_representors[qf-1];
            
            A_F(qa, qf) = inner_product(vec1, vec2);    
        }
    }
    rb->A_F_representor_norms.push_back(A_F);
    std::cout << "n = " << N-1 << " "<< rb->A_F_representor_norms[N-1] << std::endl ;

    timer.stop();
    std::cout << "  " << timer.elapsed() << " seconds" << std::endl << std::endl;
    
}

template <typename T, typename Basis, typename Prec, typename TruthSolver, typename Compression>
void
AdaptiveRBTruth2D<T, Basis, Prec, TruthSolver, Compression>::write_riesz_representors(const std::string& directory_name)
{
  // Make a directory to store all the data files
  if( mkdir(directory_name.c_str(), 0777) == -1)
  {
    std::cerr << "In RBModel::write_RB_data, directory "
                 << directory_name << " already exists, overwriting contents." << std::endl;
  }
  
  for(unsigned int i = 0; i < rb->Q_f(); ++i){
    std::stringstream filename;
    filename << directory_name << "/F_representor_" << i+1 << ".dat";
    saveCoeffVector2D(F_representors[i], basis, filename.str().c_str());
  }
  
  for(unsigned int n = 0; n < rb->n_bf(); ++n){
    for(unsigned int i = 0; i < rb->Q_a(); ++i){
      std::stringstream filename;
      filename << directory_name << "/A_representor_" << i+1 << "_" << n+1 << ".dat";
      saveCoeffVector2D(A_representors[n][i], basis, filename.str().c_str()); 
    }
  }
  
}

template <typename T, typename Basis, typename Prec, typename TruthSolver, typename Compression>
T
AdaptiveRBTruth2D<T, Basis, Prec, TruthSolver, Compression>::inner_product(const CoeffVector& v1, const CoeffVector& v2)
{
    T val = 0;
    
    if(use_inner_product_matrix && assembled_inner_product_matrix){
        // Assumption here: both vectors and the matrix have the same indexset
        
        assert(v1.size() == v2.size());
        typename CoeffVector::const_iterator it1, it2;
        
        // Build dense vectors
        DenseVectorT v1_dense(v1.size()), v2_dense(v2.size());
        int index_count = 1;
        for (it1 = v1.begin(), it2 = v2.begin(); it1 != v1.end(); ++it1, ++it2, ++index_count) {
            v1_dense(index_count) = (*it1).second;
            v2_dense(index_count) = (*it2).second;
            
        }
        
        DenseVectorT I_v2 = inner_product_matrix * v2_dense;
        val = v1_dense * I_v2;
    }
    else{
        typename CoeffVector::const_iterator it1, it2;
        for (it1 = v1.begin(); it1 != v1.end() ; ++it1) {
            for (it2 = v2.begin(); it2 != v2.end(); ++it2) {
                val += (*it1).second * (*inner_product_op)((*it1).first, (*it2).first) * (*it2).second;
            }
        }
    }
    return val;   
}

template <typename T, typename Basis, typename Prec, typename TruthSolver, typename Compression>
void
AdaptiveRBTruth2D<T, Basis, Prec, TruthSolver, Compression>::assemble_inner_product_matrix(IndexSet<Index2D>& indexset)
{
  
  Timer timer;
  std::cout << "Assemble Inner Product Matrix ...." << std::endl;
  
  inner_product_matrix.resize(indexset.size(), indexset.size());
  
  timer.start();
  
  typedef typename IndexSet<Index2D>::const_iterator const_set_it;
  std::map<Index2D,int,lt<Lexicographical,Index2D> > row_indices;
  int row_count = 1, col_count = 1;

  if(!assembled_prec_vec){
      prec_vec.engine().resize((int)indexset.size());
      for (const_set_it row=indexset.begin(); row!=indexset.end(); ++row, ++row_count) {
          row_indices[(*row)] = row_count;
          prec_vec(row_count) = get_prec((*row));
      }
      assembled_prec_vec = true;
  }
  else{
      for (const_set_it row=indexset.begin(); row!=indexset.end(); ++row, ++row_count) {
          row_indices[(*row)] = row_count;
      }    
  }
  
  repr_lhs_op.compression.setParameters(indexset);
  for (const_set_it col=indexset.begin(); col!=indexset.end(); ++col, ++col_count) {
    IndexSet<Index2D> LambdaRowSparse = repr_lhs_op.compression.SparsityPattern(*col, indexset);
    for (const_set_it row=LambdaRowSparse.begin(); row!=LambdaRowSparse.end(); ++row, ++row_count) {
        T tmp = (*inner_product_op)(*row,*col);
        if (fabs(tmp)>0)                inner_product_matrix(row_indices[*row],col_count) = tmp;
    }
  }
  inner_product_matrix.finalize();
    
  timer.stop();
  std::cout << "... done: " << timer.elapsed() << " seconds" << std::endl;
  
  assembled_inner_product_matrix = true;
}

template <typename T, typename Basis, typename Prec, typename TruthSolver, typename Compression>
void
AdaptiveRBTruth2D<T, Basis, Prec, TruthSolver, Compression>::assemble_A_operator_matrices(IndexSet<Index2D>& indexset)
{
  Timer timer;
  std::cout << "Assemble A Matrices ...." << std::endl;
  unsigned int Q_a = get_rb_model().Q_a();
  int N = (int)indexset.size();
  
  timer.start();
  for(unsigned int qa = 0; qa < Q_a; ++qa){
    SparseMatrixT A_matrix(N, N);
    
    typedef typename IndexSet<Index2D>::const_iterator const_set_it;
    std::map<Index2D,int,lt<Lexicographical,Index2D> > row_indices;
    int row_count = 1, col_count = 1;
    if((!assembled_prec_vec) && qa == 0){
        prec_vec.engine().resize((int)indexset.size());
        for (const_set_it row=indexset.begin(); row!=indexset.end(); ++row, ++row_count) {
            row_indices[(*row)] = row_count;
            prec_vec(row_count) = get_prec((*row));
        }
        assembled_prec_vec = true;
    }
    else{
        for (const_set_it row=indexset.begin(); row!=indexset.end(); ++row, ++row_count) {
            row_indices[(*row)] = row_count;
        }    
    }
    for (const_set_it col=indexset.begin(); col!=indexset.end(); ++col, ++col_count) {
        for (const_set_it row=indexset.begin(); row!=indexset.end(); ++row, ++row_count) {
            T tmp = (*A_operators[qa])(*row,*col);
            if (fabs(tmp)>0)                A_matrix(row_indices[*row],col_count) = tmp;
        }
    }
    A_matrix.finalize();

    A_operator_matrices.push_back(A_matrix);
  }
  timer.stop();
  std::cout << "... done: " << timer.elapsed() << " seconds" << std::endl;
  
  assembled_A_operator_matrices = true;
}

template <typename T, typename Basis, typename Prec, typename TruthSolver, typename Compression>
T
AdaptiveRBTruth2D<T, Basis, Prec, TruthSolver, Compression>::uncached_residual_dual_norm(const DenseVectorT& u_RB, const std::vector<T>& mu)
{
  // Calculate residual
  Operator_Residual_Representor op_res_repr(this, mu, u_RB);

  // Solve for Riesz representor
  Coefficients<Lexicographical,T,Index2D> res_repr = solver->repr_solve_totalRes(op_res_repr);
  
  // Calculate Norm
  undo_prec(res_repr);
  
  return std::sqrt(inner_product(res_repr, res_repr));
}

template <typename T, typename Basis, typename Prec, typename TruthSolver, typename Compression>
T
AdaptiveRBTruth2D<T, Basis, Prec, TruthSolver, Compression>::
uncached_residual_dual_norm(const DenseVectorT& u_RB, const std::vector<T>& mu, Coefficients<Lexicographical,T,Index2D>& res_repr)
{
  // Calculate residual
  Operator_Residual_Representor op_res_repr(this, mu, u_RB);

  // Solve for Riesz representor
  res_repr = solver->repr_solve_totalRes(op_res_repr);
  
  // Calculate Norm
  undo_prec(res_repr);
  return std::sqrt(inner_product(res_repr, res_repr));
}

// ================================================================================================================ //
// ================================================================================================================ //

/*  Operator LHS */

template <typename T, typename Basis, typename Prec, typename TruthSolver, typename Compression>
T
AdaptiveRBTruth2D<T, Basis, Prec, TruthSolver, Compression>::Operator_LHS::
operator()(const Index2D &row_index, const Index2D &col_index)
{
    if(qa < 0){
      T val = 0;
      for (unsigned int i = 0; i < thisTruth->A_operators.size(); ++i) {
          val += (*thisTruth->rb->theta_a[i])(thisTruth->rb->get_current_param()) 
               * (*thisTruth->A_operators[i])(row_index, col_index);
      }
      return thisTruth->get_prec(col_index) * thisTruth->get_prec(row_index) * val;
    }
    else{
      return thisTruth->get_prec(col_index) * thisTruth->get_prec(row_index) * (*thisTruth->A_operators[qa])(row_index, col_index);
    }
}

template <typename T, typename Basis, typename Prec, typename TruthSolver, typename Compression>
Coefficients<Lexicographical,T,Index2D>
AdaptiveRBTruth2D<T, Basis, Prec, TruthSolver, Compression>::Operator_LHS::
mv(const IndexSet<Index2D> &LambdaRow, const Coefficients<Lexicographical,T,Index2D> &x)
{
    std::cerr << "AdaptiveRBTruth2D<T, Basis, TruthSolver>::Operator_LHS::mv not implemented."
              << std::endl;
    assert(0);
    exit(1);
}

template <typename T, typename Basis, typename Prec, typename TruthSolver, typename Compression>
void
AdaptiveRBTruth2D<T, Basis, Prec, TruthSolver, Compression>::Operator_LHS::
toFlensSparseMatrix(const IndexSet<Index2D> &LambdaRow,const IndexSet<Index2D> &LambdaCol,
                    SparseMatrixT &A, T tol)
{
    std::cerr << "AdaptiveRBTruth2D<T, Basis, TruthSolver>::Operator_LHS::toFlensSparseMatrix "
              << "not implemented."
              << std::endl;
    assert(0);
    exit(1);
}

template <typename T, typename Basis, typename Prec, typename TruthSolver, typename Compression>
void
AdaptiveRBTruth2D<T, Basis, Prec, TruthSolver, Compression>::Operator_LHS::
apply(const Coefficients<Lexicographical,T,Index2D> &v, T eps,
      const IndexSet<Index2D> &Lambda, Coefficients<Lexicographical,T,Index2D> &ret,
      cxxblas::Transpose trans)
{
    std::cerr << "AdaptiveRBTruth2D<T, Basis, TruthSolver>::Operator_LHS::apply "
              << "not implemented."
              << std::endl;
    assert(0);
    exit(1);
}

/*  Operator RHS */

template <typename T, typename Basis, typename Prec, typename TruthSolver, typename Compression>
T
AdaptiveRBTruth2D<T, Basis, Prec, TruthSolver, Compression>::Operator_RHS::operator()(const Index2D &lambda)
{
    T val = 0;
    for (unsigned int i = 0; i < thisTruth->F_operators.size(); ++i) {
        val += (*thisTruth->rb->theta_f[i])(thisTruth->rb->get_current_param()) 
             * (*thisTruth->F_operators[i])(lambda);
    }
    
    return thisTruth->get_prec(lambda) * val;
}

template <typename T, typename Basis, typename Prec, typename TruthSolver, typename Compression>
Coefficients<Lexicographical,T,Index2D>
AdaptiveRBTruth2D<T, Basis, Prec, TruthSolver, Compression>::Operator_RHS::operator()(const IndexSet<Index2D> &Lambda)
{
    Coefficients<Lexicographical, T, Index2D> coeffs;
    
    typename IndexSet<Index2D>::const_iterator it_Lambda;
    typedef typename Coefficients<Lexicographical,T,Index2D>::value_type val_type;
    for (it_Lambda = Lambda.begin(); it_Lambda != Lambda.end(); ++it_Lambda) {
        coeffs.insert(val_type((*it_Lambda), this->operator()(*it_Lambda)));
    }
    
    return coeffs;
}

template <typename T, typename Basis, typename Prec, typename TruthSolver, typename Compression>
Coefficients<Lexicographical,T,Index2D>
AdaptiveRBTruth2D<T, Basis, Prec, TruthSolver, Compression>::Operator_RHS::operator()(T tol)
{
    Coefficients<Lexicographical,T,Index2D> c;
    
    std::cerr << "Operator()(T tol) : Not implemented yet" << std::endl;
    exit(1);
    
    
    return c;
}

/*  Operator LHS_Representor */

template <typename T, typename Basis, typename Prec, typename TruthSolver, typename Compression>
T
AdaptiveRBTruth2D<T, Basis, Prec, TruthSolver, Compression>::Operator_LHS_Representor::
operator()(const Index2D &row_index, const Index2D &col_index)
{
    return thisTruth->get_prec(row_index) * thisTruth->get_prec(col_index) * (*thisTruth->inner_product_op)(row_index, col_index);
}

template <typename T, typename Basis, typename Prec, typename TruthSolver, typename Compression>
Coefficients<Lexicographical,T,Index2D>
AdaptiveRBTruth2D<T, Basis, Prec, TruthSolver, Compression>::Operator_LHS_Representor::
mv(const IndexSet<Index2D> &LambdaRow, const Coefficients<Lexicographical,T,Index2D> &x)
{
    std::cerr << "AdaptiveRBTruth2D<T, Basis, TruthSolver>::Operator_LHS_Representor::mv "
              << "not implemented." << std::endl;
    assert(0);
    exit(1);
}

template <typename T, typename Basis, typename Prec, typename TruthSolver, typename Compression>
void
AdaptiveRBTruth2D<T, Basis, Prec, TruthSolver, Compression>::Operator_LHS_Representor::
toFlensSparseMatrix(const IndexSet<Index2D> &LambdaRow, const IndexSet<Index2D> &LambdaCol,
                    SparseMatrixT &A, T tol)
{
    std::cerr << "AdaptiveRBTruth2D<T, Basis, TruthSolver>::Operator_LHS_Representor::"
              <<"toFlensSparseMatrix not implemented." << std::endl;
    assert(0);
    exit(1);
}

template <typename T, typename Basis, typename Prec, typename TruthSolver, typename Compression>
void
AdaptiveRBTruth2D<T, Basis, Prec, TruthSolver, Compression>::Operator_LHS_Representor::
apply(const Coefficients<Lexicographical,T,Index2D> &v, T eps,
      const IndexSet<Index2D> &Lambda, Coefficients<Lexicographical,T,Index2D> &ret,
      cxxblas::Transpose trans)
{
    std::cerr << "AdaptiveRBTruth2D<T, Basis, TruthSolver>::Operator_LHS_Representor::"
              <<"apply not implemented." << std::endl;
    assert(0);
    exit(1);
}


/* Operator RHS_BilFormRepresentor */

template <typename T, typename Basis, typename Prec, typename TruthSolver, typename Compression>
T
AdaptiveRBTruth2D<T, Basis, Prec, TruthSolver, Compression>::Operator_RHS_BilFormRepresentor::operator()(const Index2D &lambda)
{
    T val = 0;
    typename Coefficients<Lexicographical,T,Index2D>::const_iterator it;
    for (it = current_bf->begin(); it != current_bf->end(); ++it) {
        val += (*it).second * (*current_op)((*it).first, lambda) ;
    }
    
    return - thisTruth->get_prec(lambda) * val;
}

template <typename T, typename Basis, typename Prec, typename TruthSolver, typename Compression>
Coefficients<Lexicographical,T,Index2D>
AdaptiveRBTruth2D<T, Basis, Prec, TruthSolver, Compression>::Operator_RHS_BilFormRepresentor::operator()(const IndexSet<Index2D> &Lambda)
{
    Coefficients<Lexicographical, T, Index2D> coeffs;
    
    typename IndexSet<Index2D>::const_iterator it_Lambda;
    typedef typename Coefficients<Lexicographical,T,Index2D>::value_type val_type;
    for (it_Lambda = Lambda.begin(); it_Lambda != Lambda.end(); ++it_Lambda) {
        coeffs.insert(val_type((*it_Lambda), this->operator()(*it_Lambda)));
    }
    
    return coeffs;
}

template <typename T, typename Basis, typename Prec, typename TruthSolver, typename Compression>
Coefficients<Lexicographical,T,Index2D>
AdaptiveRBTruth2D<T, Basis, Prec, TruthSolver, Compression>::Operator_RHS_BilFormRepresentor::operator()(T tol)
{
    Coefficients<Lexicographical, T, Index2D> coeffs;
    
    std::cerr << "Operator()(T tol) : Not implemented yet" << std::endl;
    exit(1);
    
    return coeffs;
}

template <typename T, typename Basis, typename Prec, typename TruthSolver, typename Compression>
void
AdaptiveRBTruth2D<T, Basis, Prec, TruthSolver, Compression>::Operator_RHS_BilFormRepresentor::
set_current_op(Operator2D<T>& op)
{
    current_op = &op;
}

template <typename T, typename Basis, typename Prec, typename TruthSolver, typename Compression>
void
AdaptiveRBTruth2D<T, Basis, Prec, TruthSolver, Compression>::Operator_RHS_BilFormRepresentor::
set_current_bf(Coefficients<Lexicographical,T,Index2D>& bf)
{
    current_bf = &bf;
}

/* Operator RHS_FunctionalRepresentor */

template <typename T, typename Basis, typename Prec, typename TruthSolver, typename Compression>
T
AdaptiveRBTruth2D<T, Basis, Prec, TruthSolver, Compression>::Operator_RHS_FunctionalRepresentor::operator()(const Index2D &lambda)
{
    return thisTruth->get_prec(lambda) * (*current_op)(lambda);
}

template <typename T, typename Basis, typename Prec, typename TruthSolver, typename Compression>
Coefficients<Lexicographical,T,Index2D>
AdaptiveRBTruth2D<T, Basis, Prec, TruthSolver, Compression>::Operator_RHS_FunctionalRepresentor::operator()(const IndexSet<Index2D> &Lambda)
{
    Coefficients<Lexicographical, T, Index2D> coeffs;
    
    typename IndexSet<Index2D>::const_iterator it_Lambda;
    typedef typename Coefficients<Lexicographical,T,Index2D>::value_type val_type;
    for (it_Lambda = Lambda.begin(); it_Lambda != Lambda.end(); ++it_Lambda) {
        coeffs.insert(val_type((*it_Lambda), this->operator()(*it_Lambda)));
    }
    
    return coeffs;
}

template <typename T, typename Basis, typename Prec, typename TruthSolver, typename Compression>
Coefficients<Lexicographical,T,Index2D>
AdaptiveRBTruth2D<T, Basis, Prec, TruthSolver, Compression>::Operator_RHS_FunctionalRepresentor::operator()(T tol)
{
    std::cerr << "Operator()(T tol) : Not implemented yet" << std::endl;
    exit(1);
    
    return (*current_op)(tol);
}

template <typename T, typename Basis, typename Prec, typename TruthSolver, typename Compression>
void
AdaptiveRBTruth2D<T, Basis, Prec, TruthSolver, Compression>::Operator_RHS_FunctionalRepresentor::
set_current_op(AdaptiveRhs<T, Index2D>& op)
{
    current_op = &op;
}

/* Operator Residual_Representor */

template <typename T, typename Basis, typename Prec, typename TruthSolver, typename Compression>
T
AdaptiveRBTruth2D<T, Basis, Prec, TruthSolver, Compression>::Operator_Residual_Representor::operator()(const Index2D &lambda)
{
  T val = 0;

  typename Coefficients<Lexicographical,T,Index2D>::const_iterator it;

  for (unsigned int i = 0; i < thisTruth->F_operators.size(); ++i) {
      val += (*thisTruth->rb->theta_f[i])(eval_mu) * (*thisTruth->F_operators[i])(lambda);
  }
   

  for (unsigned int i = 0; i < thisTruth->A_operators.size(); ++i) {
    T val_Au = 0;
    
    for (int j = 0; j < u_N.length();++j){
      T val_A_bf = 0;
      for (it = thisTruth->rb->rb_basis_functions[j].begin(); it != thisTruth->rb->rb_basis_functions[j].end(); ++it) {
          //val_A_bf += (*it).second *(*thisTruth->A_operators[i])((*it).first,lambda);
          val_A_bf += (*it).second *(*thisTruth->A_operators[i])(lambda, (*it).first);
      }
      val_Au += u_N(j+1) * val_A_bf;
    }
    
    val -= (*thisTruth->rb->theta_a[i])(eval_mu) * val_Au;
  }
  
  return thisTruth->get_prec(lambda) * val;
}

template <typename T, typename Basis, typename Prec, typename TruthSolver, typename Compression>
Coefficients<Lexicographical,T,Index2D>
AdaptiveRBTruth2D<T, Basis, Prec, TruthSolver, Compression>::Operator_Residual_Representor::operator()(const IndexSet<Index2D> &Lambda)
{
  Coefficients<Lexicographical, T, Index2D> coeffs;
  
  typename IndexSet<Index2D>::const_iterator it_Lambda;
  typedef typename Coefficients<Lexicographical,T,Index2D>::value_type val_type;
  for (it_Lambda = Lambda.begin(); it_Lambda != Lambda.end(); ++it_Lambda) {
      coeffs.insert(val_type((*it_Lambda), this->operator()(*it_Lambda)));
  }
  
  return coeffs;
}

template <typename T, typename Basis, typename Prec, typename TruthSolver, typename Compression>
Coefficients<Lexicographical,T,Index2D>
AdaptiveRBTruth2D<T, Basis, Prec, TruthSolver, Compression>::Operator_Residual_Representor::operator()(T tol)
{
  Coefficients<Lexicographical, T, Index2D> coeffs;
  
  std::cerr << "Operator()(T tol) : Not implemented yet" << std::endl;
  exit(1);
  
  return coeffs;
}



} // namespace lawa
