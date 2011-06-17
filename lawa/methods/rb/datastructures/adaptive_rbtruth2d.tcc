#include <lawa/methods/rb/postprocessing/plotting.h>

namespace  lawa {

template <typename T, typename Basis, typename TruthSolver, typename Compression>
AdaptiveRBTruth2D<T, Basis, TruthSolver, Compression>::AdaptiveRBTruth2D(Basis& _basis, bool _use_inner_product, 
                      bool _use_A_matrix, bool _use_F_vector)
    : basis(_basis), lhs_op(this), rhs_op(this), repr_lhs_op(this), repr_rhs_A_op(this), repr_rhs_F_op(this),
      use_inner_product_matrix(_use_inner_product), use_A_operator_matrices(_use_A_matrix), 
      use_F_operator_vectors(_use_F_vector)
{}

template <typename T, typename Basis, typename TruthSolver, typename Compression>
void
AdaptiveRBTruth2D<T, Basis, TruthSolver, Compression>::attach_A_q(theta_fctptr theta_a_q, Operator2D<T>& A_q)
{
    rb->theta_a.push_back(theta_a_q);
    A_operators.push_back(&A_q);
    
    // If we have a fixed index set, compute the corresponding matrix
   /* if(flens::IsSame<IndexsetSolver, TruthSolver>::value){
      Timer timer;
      SparseMatrixT A_q_matrix;
      std::cout << "Compute Operator matrix .... " << std::endl;
      timer.start();
      toFlensSparseMatrix(A_q, solver->basis_set, solver->basis_set, A_q_matrix);
      timer.stop();
      A_operator_matrices.push_back(A_q_matrix);
      std::cout << ".... done: " << timer.elapsed() << " seconds" << std::endl;
    }*/
}

template <typename T, typename Basis, typename TruthSolver, typename Compression>
void
AdaptiveRBTruth2D<T, Basis, TruthSolver, Compression>::attach_A_q(Operator2D<T>& A_q)
{
    A_operators.push_back(&A_q);
}

template <typename T, typename Basis, typename TruthSolver, typename Compression>
void
AdaptiveRBTruth2D<T, Basis, TruthSolver, Compression>::attach_F_q(theta_fctptr theta_f_q, AdaptiveRhs<T, Index2D>& F_q)
{
    rb->theta_f.push_back(theta_f_q);
    F_operators.push_back(&F_q);
    
    // If we have a fixed index set, compute the corresponding vector
   /* if(flens::IsSame<IndexsetSolver, TruthSolver>::value){
      Timer timer;
      DenseVectorT F_q_vector(solver->basis_set.size());
      
      typedef typename IndexSet<Index >::const_iterator const_set_it;
      
      std::cout << "Compute RHS vector .... " << std::endl;
      timer.start();
      int count=1;
      for (const_set_it it = solver->basis_set.begin(); it != solver->basis_set.end(); ++it, ++count) {
        F_q_vector(count) = F_q((*it));
      }
      timer.stop();
      F_operator_vectors.push_back(F_q_vector);
      std::cout << ".... done: " << timer.elapsed() << " seconds" << std::endl;
    }*/
}

template <typename T, typename Basis, typename TruthSolver, typename Compression>
void
AdaptiveRBTruth2D<T, Basis, TruthSolver, Compression>::attach_F_q(AdaptiveRhs<T, Index2D>& F_q)
{
    F_operators.push_back(&F_q);
}

template <typename T, typename Basis, typename TruthSolver, typename Compression>
void 
AdaptiveRBTruth2D<T, Basis, TruthSolver, Compression>::set_truthsolver(TruthSolver& _truthsolver)
{
    solver = &_truthsolver;
    solver->set_model(*this);
    
    // If we have a fixed index set, compute the inner product matrix
    /*if(flens::IsSame<IndexsetSolver, TruthSolver>::value){
      Timer timer;
      std::cout << "Compute Inner product matrix .... " << std::endl;
      timer.start();
      toFlensSparseMatrix(*rb->inner_product_op, solver->basis_set, solver->basis_set, rb->inner_product_matrix);
      timer.stop();
      std::cout << ".... done: " << timer.elapsed() << " seconds" << std::endl;
    }*/
}

template <typename T, typename Basis, typename TruthSolver, typename Compression>
void
AdaptiveRBTruth2D<T, Basis, TruthSolver, Compression>::set_rb_model(RBModel2D<T, AdaptiveRBTruth2D<T, Basis, TruthSolver, Compression> >& _rb)
{
    rb = &_rb;
}

template <typename T, typename Basis, typename TruthSolver, typename Compression>
RBModel2D<T, AdaptiveRBTruth2D<T, Basis, TruthSolver, Compression> >&
AdaptiveRBTruth2D<T, Basis, TruthSolver, Compression>::get_rb_model()
{
    return *rb;
}

template <typename T, typename Basis, typename TruthSolver, typename Compression>
void
AdaptiveRBTruth2D<T, Basis, TruthSolver, Compression>::update_representors()
{
    //Timer timer;
    int N = rb->n_bf();
    if(N == 1) {
        F_representors.resize(rb->Q_f());
        for (unsigned int i = 0; i < rb->Q_f(); ++i) {
            repr_rhs_F_op.set_current_op(*F_operators[i]);
            std::cout<< std::endl << "  --- Solving for Riesz Representor of F_" << i+1 << "---" << std::endl<< std::endl;
            CoeffVector c = solver->repr_solve_F();
            
            /*
            // Scatterplot
            std::cout << "Plotting Representor F_"<< i+1 << " .... " << std::endl;
            std::stringstream s;
            s << "F_Representor_" << i+1;
            timer.start();
            saveCoeffVector2D(c, basis, s.str().c_str());
            timer.stop();
            std::cout << ".... finished: " << timer.elapsed() << " seconds" << std::endl << std::endl;
            */
            
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
        
       /*
        // Scatterplot
        std::cout << "  Plotting Representor A_"<< i+1 << " .... " << std::endl;
        std::stringstream s;
        s << "A_Representor_" << i+1 << "_" << N;
        timer.start();
        saveCoeffVector2D(c, basis, s.str().c_str());
        timer.stop();
        std::cout << "  .... finished: " << timer.elapsed() << " seconds" << std::endl << std::endl;
       */
            
        A_representors[N-1][i] = c;
    }
}

template <typename T, typename Basis, typename TruthSolver, typename Compression>
void
AdaptiveRBTruth2D<T, Basis, TruthSolver, Compression>::update_representor_norms()
{
  Timer timer;
   // Coefficients<Lexicographical,T,Index> X_F_or_A;
    int N = rb->n_bf();
    int Qf = rb->Q_f();

    if(N == 1){
      std::cout << "  ... F x F .... " << std::endl;
       // F_F representor norms 
       rb->F_F_representor_norms.engine().resize(Qf, Qf);

       timer.start();
       for(int qf1 = 1; qf1 <= Qf; ++qf1) {
           for (int qf2 = qf1; qf2 <= Qf; ++qf2) {
               //std::cout << ".... F_" << qf1 << " ,  F_" << qf2 << std::endl;
               //X_F_or_A = mv_sparse(supp(F_representors[qf2-1]),rb->inner_product,F_representors[qf2-1]);
               //F_F_representor_norms(qf1,qf2) = F_representors[qf1-1]*X_F_or_A;

               rb->F_F_representor_norms(qf1,qf2) = rb->inner_product(F_representors[qf1-1], F_representors[qf2-1]);

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
        //  n2 = N
        flens::GeMatrix<flens::FullStorage<T, cxxblas::ColMajor> > A_n1_N(Qa, Qa);
        for(int qa1 = 1; qa1 <= Qa; ++qa1) {
          for(int qa2 = qa1; qa2 <= Qa; ++qa2) {                  
            A_n1_N(qa1, qa2) = rb->inner_product(A_representors[n1][qa1-1], A_representors[N-1][qa2-1]);
            if(qa1 != qa2){
              if(n1 == N-1){
                A_n1_N(qa2, qa1) = A_n1_N(qa1, qa2);
              }
              else{
                A_n1_N(qa2, qa1) = rb->inner_product(A_representors[n1][qa2-1], A_representors[N-1][qa1-1]);                              
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

    std::cout << "  ... A x F .... " << std::endl;
    // A_F representor norms 
    timer.start();
    // n = N
    flens::GeMatrix<flens::FullStorage<T, cxxblas::ColMajor> > A_F(Qa, Qf);
    for(int qa = 1; qa <= Qa; ++qa) {
        for(int qf = 1; qf <= Qf; ++qf) {
           // std::cout << "....A_" << qa <<"(" << n << "), F_" << qf << std::endl;
            A_F(qa, qf) = rb->inner_product(A_representors[N-1][qa-1], F_representors[qf-1]);    
        }
    }
    rb->A_F_representor_norms.push_back(A_F);
    std::cout << "n = " << N-1 << " "<< rb->A_F_representor_norms[N-1] << std::endl ;

    timer.stop();
    std::cout << "  " << timer.elapsed() << " seconds" << std::endl << std::endl;
    
}

template <typename T, typename Basis, typename TruthSolver, typename Compression>
void
AdaptiveRBTruth2D<T, Basis, TruthSolver, Compression>::write_riesz_representors(const std::string& directory_name)
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

template <typename T, typename Basis, typename TruthSolver, typename Compression>
void
AdaptiveRBTruth2D<T, Basis, TruthSolver, Compression>::assemble_inner_product_matrix(IndexSet<Index2D>& indexset)
{
  Timer timer;
  std::cout << "Assemble Inner Product Matrix ...." << std::endl;
  inner_product_matrix.resize(indexset.size(), indexset.size());
  timer.start();
  toFlensSparseMatrix(repr_lhs_op, indexset, indexset, inner_product_matrix);
  timer.stop();
  std::cout << "... done: " << timer.elapsed() << " seconds" << std::endl;

  rb->assembled_inner_product_matrix = true;
}

template <typename T, typename Basis, typename TruthSolver, typename Compression>
void
AdaptiveRBTruth2D<T, Basis, TruthSolver, Compression>::assemble_A_operator_matrices(IndexSet<Index2D>& indexset)
{
  Timer timer;
  std::cout << "Assemble A Matrices ...." << std::endl;
  unsigned int Q_a = get_rb_model().Q_a();
  int N = indexset.size();
  
  timer.start();
  for(unsigned int qa = 1; qa <= Q_a; ++qa){
    SparseMatrixT A_matrix(N, N);
    lhs_op.qa = qa-1;
    toFlensSparseMatrix(lhs_op, indexset, indexset, A_matrix);
    A_operator_matrices.push_back(A_matrix);
  }
  lhs_op.qa = -1;
  timer.stop();
  std::cout << "... done: " << timer.elapsed() << " seconds" << std::endl;
  
  rb->assembled_A_operator_matrices = true;
}
// ================================================================================================================ //
// ================================================================================================================ //

/*  Operator LHS */

template <typename T, typename Basis, typename TruthSolver, typename Compression>
T
AdaptiveRBTruth2D<T, Basis, TruthSolver, Compression>::Operator_LHS::
operator()(const Index2D &row_index, const Index2D &col_index)
{
    if(qa < 0){
      T val = 0;
      for (unsigned int i = 0; i < thisTruth->A_operators.size(); ++i) {
          val += (*thisTruth->rb->theta_a[i])(thisTruth->rb->get_current_param()) 
               * (*thisTruth->A_operators[i])(row_index, col_index);
      }
      return val;
    }
    else{
      return (*thisTruth->A_operators[qa])(row_index, col_index);
    }
}

/*  Operator RHS */

template <typename T, typename Basis, typename TruthSolver, typename Compression>
T
AdaptiveRBTruth2D<T, Basis, TruthSolver, Compression>::Operator_RHS::operator()(const Index2D &lambda)
{
    T val = 0;
    for (unsigned int i = 0; i < thisTruth->F_operators.size(); ++i) {
        val += (*thisTruth->rb->theta_f[i])(thisTruth->rb->get_current_param()) 
             * (*thisTruth->F_operators[i])(lambda);
    }
    
    return val;
}

template <typename T, typename Basis, typename TruthSolver, typename Compression>
Coefficients<Lexicographical,T,Index2D>
AdaptiveRBTruth2D<T, Basis, TruthSolver, Compression>::Operator_RHS::operator()(const IndexSet<Index2D> &Lambda)
{
    Coefficients<Lexicographical,T,Index2D> c;
    for (unsigned int i = 0; i < thisTruth->F_operators.size(); ++i) {
        c = c + (*thisTruth->F_operators[i])(Lambda) 
              * (*thisTruth->rb->theta_f[i])(thisTruth->rb->get_current_param());
    }
    
    return c;
}

template <typename T, typename Basis, typename TruthSolver, typename Compression>
Coefficients<Lexicographical,T,Index2D>
AdaptiveRBTruth2D<T, Basis, TruthSolver, Compression>::Operator_RHS::operator()(T tol)
{
    Coefficients<Lexicographical,T,Index2D> c;
    for (unsigned int i = 0; i < thisTruth->F_operators.size(); ++i) {
        c += (*thisTruth->F_operators[i])(tol) 
           * (*thisTruth->rb->theta_f[i])(thisTruth->rb->get_current_param());
    }
    
    return c;
}

/*  Operator LHS_Representor */

template <typename T, typename Basis, typename TruthSolver, typename Compression>
T
AdaptiveRBTruth2D<T, Basis, TruthSolver, Compression>::Operator_LHS_Representor::
operator()(const Index2D &row_index, const Index2D &col_index)
{
    return (*thisTruth->rb->inner_product_op)(row_index, col_index);
}

/* Operator RHS_BilFormRepresentor */

template <typename T, typename Basis, typename TruthSolver, typename Compression>
T
AdaptiveRBTruth2D<T, Basis, TruthSolver, Compression>::Operator_RHS_BilFormRepresentor::operator()(const Index2D &lambda)
{
    T val = 0;
    typename Coefficients<Lexicographical,T,Index2D>::const_iterator it;
    for (it = current_bf->begin(); it != current_bf->end(); ++it) {
        val += (*it).second * (*current_op)((*it).first, lambda);
    }
    
    return - val;
}

template <typename T, typename Basis, typename TruthSolver, typename Compression>
Coefficients<Lexicographical,T,Index2D>
AdaptiveRBTruth2D<T, Basis, TruthSolver, Compression>::Operator_RHS_BilFormRepresentor::operator()(const IndexSet<Index2D> &Lambda)
{
    Coefficients<Lexicographical, T, Index2D> coeffs;
    
    typename IndexSet<Index2D>::const_iterator it_Lambda;
    typedef typename Coefficients<Lexicographical,T,Index2D>::value_type val_type;
    for (it_Lambda = Lambda.begin(); it_Lambda != Lambda.end(); ++it_Lambda) {
        coeffs.insert(val_type((*it_Lambda), this->operator()(*it_Lambda)));
    }
    
    return coeffs;
}

template <typename T, typename Basis, typename TruthSolver, typename Compression>
Coefficients<Lexicographical,T,Index2D>
AdaptiveRBTruth2D<T, Basis, TruthSolver, Compression>::Operator_RHS_BilFormRepresentor::operator()(T tol)
{
    Coefficients<Lexicographical, T, Index2D> coeffs;
    
    std::cerr << "Operator()(T tol) : Not implemented yet" << std::endl;
    exit(1);
    
    return coeffs;
}

template <typename T, typename Basis, typename TruthSolver, typename Compression>
void
AdaptiveRBTruth2D<T, Basis, TruthSolver, Compression>::Operator_RHS_BilFormRepresentor::
set_current_op(Operator2D<T>& op)
{
    current_op = &op;
}

template <typename T, typename Basis, typename TruthSolver, typename Compression>
void
AdaptiveRBTruth2D<T, Basis, TruthSolver, Compression>::Operator_RHS_BilFormRepresentor::
set_current_bf(Coefficients<Lexicographical,T,Index2D>& bf)
{
    current_bf = &bf;
}

/* Operator RHS_FunctionalRepresentor */

template <typename T, typename Basis, typename TruthSolver, typename Compression>
T
AdaptiveRBTruth2D<T, Basis, TruthSolver, Compression>::Operator_RHS_FunctionalRepresentor::operator()(const Index2D &lambda)
{
    return (*current_op)(lambda);
}

template <typename T, typename Basis, typename TruthSolver, typename Compression>
Coefficients<Lexicographical,T,Index2D>
AdaptiveRBTruth2D<T, Basis, TruthSolver, Compression>::Operator_RHS_FunctionalRepresentor::operator()(const IndexSet<Index2D> &Lambda)
{
    return (*current_op)(Lambda);
}

template <typename T, typename Basis, typename TruthSolver, typename Compression>
Coefficients<Lexicographical,T,Index2D>
AdaptiveRBTruth2D<T, Basis, TruthSolver, Compression>::Operator_RHS_FunctionalRepresentor::operator()(T tol)
{
    return (*current_op)(tol);
}

template <typename T, typename Basis, typename TruthSolver, typename Compression>
void
AdaptiveRBTruth2D<T, Basis, TruthSolver, Compression>::Operator_RHS_FunctionalRepresentor::
set_current_op(AdaptiveRhs<T, Index2D>& op)
{
    current_op = &op;
}

} // namespace lawa
