namespace lawa {

template <typename T, typename Basis, typename Index, typename Compression>
IndexsetTruthSolver<T, Basis, Index, Compression>::
IndexsetTruthSolver(IndexSet<Index>& _indexset, Truth& _truth, SolverCall solmethod, bool _use_inner_product, 
                    bool _use_A_matrix, bool _use_F_vector, T _tol, int maxIts)
  : basis_set(_indexset), truth_model(&_truth), solution_method(solmethod), 
    use_inner_product_matrix(_use_inner_product), use_A_operator_matrices(_use_A_matrix), 
    use_F_operator_vectors(_use_F_vector), tol(_tol), maxIterations(maxIts)
{}

template <typename T, typename Basis, typename Index, typename Compression>
void 
IndexsetTruthSolver<T, Basis, Index, Compression>
::set_model(AdaptiveRBTruth2D<T, Basis, IndexsetTruthSolver<T, Basis, Index, Compression>, Compression >& _truth_model){
    truth_model = &_truth_model;
}

template <typename T, typename Basis, typename Index, typename Compression>
Coefficients<Lexicographical,T,Index>
IndexsetTruthSolver<T, Basis, Index, Compression>::truth_solve()
{ 
  Coefficients<Lexicographical,T,Index> u, f;
  std::cout << "Start Truth Solve: Call Method " << solution_method << std::endl; 
  
  Timer timer;
  // Construct rhs coefficients vector
  typename IndexSet<Index>::const_iterator it;
  typedef typename Coefficients<Lexicographical,T,Index2D>::value_type val_type;
  std::cout << "  Build Rhs... " <<  std::endl;
  timer.start(); 
  for(it = basis_set.begin(); it != basis_set.end(); ++it){
    f.insert(val_type((*it), truth_model->rhs_op(*it)));
  }
  timer.stop();
  std::cout << "  ..... done : " << timer.elapsed() << " seconds" << std::endl;
  
  T res;
  int its;
  timer.start();
  switch(solution_method){
    case call_cg:
    std::cout << "  Start CG Solve: Maximal iterations = " << maxIterations << std::endl; 
      its = CG_Solve(basis_set, truth_model->lhs_op, u, f, res, tol, maxIterations);
      std::cout << "  CG iterations: " << its << ", residual = " << res << std::endl;
      break;
    case call_gmres:
      its = GMRES_Solve(basis_set, truth_model->lhs_op, u, f, res, tol, maxIterations);
      std::cout << "  GMRES iterations: " << its << ", residual = " << res << std::endl;
      break;
    default: 
      std::cerr << "Method not implemented yet " << std::endl;
      break;
  }
  timer.stop();
  std::cout << "Done Truth Solve: "<< timer.elapsed() << " seconds" << std::endl << std::endl;
  
  
  return u;
}

template <typename T, typename Basis, typename Index, typename Compression>
Coefficients<Lexicographical,T,Index>
IndexsetTruthSolver<T, Basis, Index, Compression>::repr_solve_F()
{
  Coefficients<Lexicographical,T,Index> u, f;
  Timer timer;
  
  if(!use_inner_product_matrix){ // Assemble LHS Matrix
  
    // Construct rhs coefficients vector
    typename IndexSet<Index>::const_iterator it;
    typedef typename Coefficients<Lexicographical,T,Index2D>::value_type val_type;
    std::cout << "  Build Rhs... " <<  std::endl;
    timer.start();
    for(it = basis_set.begin(); it != basis_set.end(); ++it){
      f.insert(val_type((*it), truth_model->repr_rhs_F_op(*it)));
    }
    timer.stop();
    std::cout << "  .... done : " << timer.elapsed() << " seconds" << std::endl;
  
    T res;
    int its;
    timer.start();
    switch(solution_method){
      case call_cg:
      std::cout << "  Start CG Solve: Maximal iterations = " << maxIterations << std::endl; 
        its = CG_Solve(basis_set, truth_model->repr_lhs_op, u, f, res, tol, maxIterations);
        std::cout << "  CG iterations: " << its << ", residual = " << res << std::endl;
        break;
      case call_gmres:
        its = GMRES_Solve(basis_set, truth_model->repr_lhs_op, u, f, res, tol, maxIterations);
        std::cout << "  GMRES iterations: " << its << ", residual = " << res << std::endl;
        break;
      default: 
        std::cerr << "Method not implemented yet " << std::endl;
        break;
    }
    timer.stop();
    std::cout << "Done F Representor Solve: "<< timer.elapsed() << " seconds" << std::endl << std::endl;
  }
  else{ // Use pre-assembled matrix (assumes assemble_inner_product_matrix has been called before)
    
    std::cout << "    Using inner product matrix!" << std::endl;
    if (basis_set.size() > 0) {
      
      typedef typename IndexSet<Index>::const_iterator const_set_it;
      typedef typename Coefficients<Lexicographical,T,Index >::const_iterator const_coeff_it;
          
      std::cout << "    Build Dense Vectors ..." << std::endl;
      
      timer.start();
      int N = basis_set.size();
      DenseVector<Array<T> > rhs(N), x(N), res(N), Ax(N);
      int row_count=1;
      
      for (const_set_it row = basis_set.begin(); row != basis_set.end(); ++row, ++row_count) {
        rhs(row_count) = truth_model->repr_rhs_F_op((*row));
        x(row_count) = 0.;  
      }
      timer.stop();
      std::cout << "    .... done : " << timer.elapsed() << " seconds " << std::endl;  
      
      T residual;
      int its;
      timer.start();
      switch(solution_method){
        case call_cg:
          std::cout << "  Start CG Solve: Maximal iterations = " << maxIterations << std::endl; 
          its = lawa::cg(inner_product_matrix, x, rhs, tol, maxIterations);
          Ax = inner_product_matrix*x;
          res= Ax-rhs;
          residual = std::sqrt(res*res);
          row_count = 1;
          for (const_set_it row=basis_set.begin(); row!=basis_set.end(); ++row, ++row_count) {
              u[*row] = x(row_count);
          }
          std::cout << "  CG iterations: " << its << ", residual = " << residual << std::endl;
          break;
        case call_gmres:
          std::cerr << "Method not implemented yet " << std::endl;
          break;
        default: 
          std::cerr << "Method not implemented yet " << std::endl;
          break;
      }
      timer.stop();
      std::cout << "Done F Representor Solve: "<< timer.elapsed() << " seconds" << std::endl << std::endl;
      
    }
  }
  
  return u;
}

template <typename T, typename Basis, typename Index, typename Compression>
Coefficients<Lexicographical,T,Index>
IndexsetTruthSolver<T, Basis, Index, Compression>::repr_solve_A()
{
  Timer timer;
  Coefficients<Lexicographical,T,Index> u;
  
  if(!use_inner_product_matrix){ // Assemble LHS Matrix
    
    Coefficients<Lexicographical,T,Index> f;  
    // Construct rhs coefficients vector
    typename IndexSet<Index>::const_iterator it;
    typedef typename Coefficients<Lexicographical,T,Index2D>::value_type val_type;
    std::cout << "  Build Rhs... " <<  std::endl;
    timer.start();
    for(it = basis_set.begin(); it != basis_set.end(); ++it){
      f.insert(val_type((*it), truth_model->repr_rhs_A_op(*it)));
    }
    timer.stop();
    std::cout << "  .... done : " << timer.elapsed() << " seconds" << std::endl;

    T res;
    int its;
    timer.start();
    switch(solution_method){
      case call_cg:
        std::cout << "  Start CG Solve: Maximal iterations = " << maxIterations << std::endl; 
        its = CG_Solve(basis_set, truth_model->repr_lhs_op, u, f, res, tol, maxIterations);
        std::cout << "  CG iterations: " << its << ", residual = " << res << std::endl;
        break;
      case call_gmres:
        its = GMRES_Solve(basis_set, truth_model->repr_lhs_op, u, f, res, tol, maxIterations);
        std::cout << "  GMRES iterations: " << its << ", residual = " << res << std::endl;
        break;
      default: 
        std::cerr << "Method not implemented yet " << std::endl;
        break;
    }
    timer.stop();
    std::cout << "Done A Representor Solve: "<< timer.elapsed() << " seconds" << std::endl << std::endl;    
  }
  else{ // Use pre-assembled matrix (assumes assemble_inner_product_matrix has been called before)
    
    std::cout << "    Using inner product matrix!" << std::endl;
    if (basis_set.size() > 0) {
      
      typedef typename IndexSet<Index>::const_iterator const_set_it;
      typedef typename Coefficients<Lexicographical,T,Index >::const_iterator const_coeff_it;
          
      std::cout << "    Build Dense Vectors ..." << std::endl;
      
      timer.start();
      int N = basis_set.size();
      DenseVector<Array<T> > rhs(N), x(N), res(N), Ax(N);
      int row_count=1;
      
      for (const_set_it row = basis_set.begin(); row != basis_set.end(); ++row, ++row_count) {
        rhs(row_count) = truth_model->repr_rhs_A_op((*row));
        x(row_count) = 0.;  
      }
      timer.stop();
      std::cout << "    .... done : " << timer.elapsed() << " seconds " << std::endl;  
      
      T residual;
      int its;
      timer.start();
      switch(solution_method){
        case call_cg:
          std::cout << "  Start CG Solve: Maximal iterations = " << maxIterations << std::endl; 
          its = lawa::cg(inner_product_matrix, x, rhs, tol, maxIterations);
          Ax = inner_product_matrix*x;
          res= Ax-rhs;
          residual = std::sqrt(res*res);
          row_count = 1;
          for (const_set_it row=basis_set.begin(); row!=basis_set.end(); ++row, ++row_count) {
              u[*row] = x(row_count);
          }
          std::cout << "  CG iterations: " << its << ", residual = " << residual << std::endl;
          break;
        case call_gmres:
          std::cerr << "Method not implemented yet " << std::endl;
          break;
        default: 
          std::cerr << "Method not implemented yet " << std::endl;
          break;
      }
      timer.stop();
      std::cout << "Done A Representor Solve: "<< timer.elapsed() << " seconds" << std::endl << std::endl;
      
    }
  }

  
  return u;  
}

template <typename T, typename Basis, typename Index, typename Compression>
void
IndexsetTruthSolver<T, Basis, Index, Compression>::assemble_inner_product_matrix()
{
  Timer timer;
  std::cout << "Assemble Inner Product Matrix ...." << std::endl;
  inner_product_matrix.resize(basis_set.size(), basis_set.size());
  timer.start();
  toFlensSparseMatrix(truth_model->repr_lhs_op, basis_set, basis_set, inner_product_matrix);
  timer.stop();
  std::cout << "... done: " << timer.elapsed() << " seconds" << std::endl;
}

template <typename T, typename Basis, typename Index, typename Compression>
void
IndexsetTruthSolver<T, Basis, Index, Compression>::assemble_A_operator_matrices()
{
  Timer timer;
  std::cout << "Assemble A Matrices ...." << std::endl;
  unsigned int Q_a = truth_model->get_rb_model().Q_a();
  int N = basis_set.size();
  
  timer.start();
  for(unsigned int qa = 1; qa <= Q_a; ++qa){
    SparseMatrixT A_matrix(N, N);
    truth_model->lhs_op.qa = qa-1;
    toFlensSparseMatrix(truth_model->lhs_op, basis_set, basis_set, A_matrix);
    A_operator_matrices.push_back(A_matrix);
  }
  truth_model->lhs_op.qa = -1;
  timer.stop();
  std::cout << "... done: " << timer.elapsed() << " seconds" << std::endl;
}

/*template <typename T, typename Basis, typename Index, typename Compression>
void
IndexsetTruthSolver<T, Basis, Index, Compression>::assemble_F_operator_vectors()
{
  Timer timer;
  std::cout << "Assemble F vectors ...." << std::endl;
  unsigned int Q_a = truth_model->get_rb_model.Q_a();
  int N = basis_set.size();
  
  timer.start();
  for(unsigned int qa = 1; qa <= Q_a; ++qa){
    DenseVectorT A_matrix(N, N);
    toFlensSparseMatrix(truth_model->A_operators[i], basis_set, basis_set, A_matrix);
    A_operator_matrices.push_back(A_matrix);
  }
  timer.stop();
  std::cout << "... done: " << timer.elapsed() << " seconds" << std::endl;
}*/

} // namespace lawa
