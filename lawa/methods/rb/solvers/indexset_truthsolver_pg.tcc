namespace lawa {

template <typename T, typename TrialBasis, typename Index, typename Compression, typename TestBasis>
IndexsetTruthSolver_PG<T, TrialBasis, Index, Compression, TestBasis>::
IndexsetTruthSolver_PG(IndexSet<Index>& _indexset_trial, IndexSet<Index>& _indexset_test, Truth& _truth, SolverCall solmethod, T _tol, int maxIts)
  : trialbasis_set(_indexset_trial), testbasis_set(_indexset_test) truth_model(&_truth), solution_method(solmethod), tol(_tol), maxIterations(maxIts)
{}

template <typename T, typename TrialBasis, typename Index, typename Compression, typename TestBasis>
void 
IndexsetTruthSolver_PG<T, TrialBasis, Index, Compression, TestBasis>
::set_model(AdaptiveRBTruth2D<T, TrialBasis, IndexsetTruthSolver_PG<T, TrialBasis, Index, Compression, TestBasis>, Compression,TrialBasis >& _truth_model){
    truth_model = &_truth_model;
}

template <typename T, typename TrialBasis, typename Index, typename Compression, typename TestBasis>
Coefficients<Lexicographical,T,Index>
IndexsetTruthSolver_PG<T, TrialBasis, Index, Compression, TestBasis>::truth_solve()
{ 
  Coefficients<Lexicographical,T,Index> u, f;
  std::cout << "Start Truth Solve: Call Method " << solution_method << std::endl; 
  
  Timer timer1, timer2;
  
  timer1.start();
  
  // Construct rhs coefficients vector
  typename IndexSet<Index>::const_iterator it;
  typedef typename Coefficients<Lexicographical,T,Index2D>::value_type val_type;
  std::cout << "  Build Rhs... " <<  std::endl;
  timer2.start(); 
  for(it = testbasis_set.begin(); it != testbasis_set.end(); ++it){
    f.insert(val_type((*it), truth_model->rhs_op(*it)));
  }
  timer2.stop();
  std::cout << "  ..... done : " << timer2.elapsed() << " seconds" << std::endl;
  
  T res;
  int its;
  switch(solution_method){
    case call_cg:
    std::cout << " This only makes sense if TrialBasis == TestBasis !! Are you sure?" << std::endl;
    if(testbasis_set.size() != trialbasis_set.size()){
      std::cerr << "Dimension of Trial and Test Basis are different -> Cannot apply CG! " << std::endl;
      exit(1);
    }
    std::cout << "  Start CG Solve: Maximal iterations = " << maxIterations << std::endl; 
      its = CG_Solve(trialbasis_set, truth_model->lhs_op, u, f, res, tol, maxIterations);
      std::cout << "  CG iterations: " << its << ", residual = " << res << std::endl;
      break;
    case call_gmres:
    if(testbasis_set.size() != trialbasis_set.size()){
      std::cerr << "Dimension of Trial and Test Basis are different -> Cannot apply GMRES! " << std::endl;
      exit(1);
    }
    std::cout << "  Start GMRES Solve: Maximal iterations = " << maxIterations << std::endl; 
      its = GMRES_Solve_PG(testbasis_set, trialbasis_set, truth_model->lhs_op, u, f, res, tol, maxIterations);
      std::cout << "  GMRES iterations: " << its << ", residual = " << res << std::endl;
      break;
    case call_cgls: 
      if(testbasis_set.size() < trialbasis_set.siz()){
        std::cerr << "Dimensions of TestBasis smaller than that of TrialBasis -> Cannot apply CGLS! " << std::endl;
        exit(1);
      }
      std::cout << "  Start CGLS Solve: Maximal iterations = " << maxIterations << std::endl; 
      its = CGLS_Solve(testbasis_set, trialbasis_set, truth_model->lhs_op, u, f, res, tol, maxIterations);
      std::cout << "  CGLS iterations: " << its << ", residual = " << res << std::endl;
      break;      
    default: 
      std::cerr << "Method not implemented yet " << std::endl;
      break;
  }
  timer1.stop();
  std::cout << "Done Truth Solve: "<< timer1.elapsed() << " seconds" << std::endl << std::endl;
  
  return u;
}

// Representor Solves: Only in Test Space 
template <typename T, typename TrialBasis, typename Index, typename Compression, typename TestBasis>
Coefficients<Lexicographical,T,Index>
IndexsetTruthSolver_PG<T, TrialBasis, Index, Compression, TestBasis>::repr_solve_F()
{
  Coefficients<Lexicographical,T,Index> u, f;
  Timer timer1, timer2;
  
  timer1.start();
  
  if(!truth_model->use_inner_product_matrix){ // Assemble LHS Matrix
  
    // Construct rhs coefficients vector
    typename IndexSet<Index>::const_iterator it;
    typedef typename Coefficients<Lexicographical,T,Index2D>::value_type val_type;
    std::cout << "  Build Rhs... " <<  std::endl;
    timer2.start();
    for(it = testbasis_set.begin(); it != testbasis_set.end(); ++it){
      f.insert(val_type((*it), truth_model->repr_rhs_F_op(*it)));
    }
    timer2.stop();
    std::cout << "  .... done : " << timer2.elapsed() << " seconds" << std::endl;
  
    T res;
    int its;
    switch(solution_method){
      case call_cg:
        std::cout << "  Start CG Solve: Maximal iterations = " << maxIterations << std::endl; 
        its = CG_Solve(testbasis_set, truth_model->repr_lhs_op, u, f, res, tol, maxIterations);
        std::cout << "  CG iterations: " << its << ", residual = " << res << std::endl;
        break;
      case call_gmres:
        its = GMRES_Solve_PG(testbasis_set, testbasis_set, truth_model->repr_lhs_op, u, f, res, tol, maxIterations);
        std::cout << "  GMRES iterations: " << its << ", residual = " << res << std::endl;
        break;
      case call_cgls:
        std::cout << "  Start CGLS Solve: Maximal iterations = " << maxIterations << std::endl; 
        its = CGLS_Solve(testbasis_set, testbasis_set, truth_model->repr_lhs_op, u, f, res, tol, maxIterations);
        std::cout << "  CGLS iterations: " << its << ", residual = " << res << std::endl;
        break;
      default: 
        std::cerr << "Method not implemented yet " << std::endl;
        break;
    }
    timer1.stop();
    std::cout << "Done F Representor Solve: "<< timer1.elapsed() << " seconds" << std::endl << std::endl;
  }
  else{ // Use pre-assembled matrix (assumes assemble_inner_product_matrix has been called before)
    
    std::cout << "    Using inner product matrix!" << std::endl;
    if (trialbasis_set.size() > 0) {
      
      typedef typename IndexSet<Index>::const_iterator const_set_it;
      typedef typename Coefficients<Lexicographical,T,Index >::const_iterator const_coeff_it;
          
      std::cout << "    Build Dense Vectors ..." << std::endl;
      
      timer2.start();
      int numRows = testbasis_set.size();
      int numCols = testbasis_set.size();
      DenseVector<Array<T> > rhs(numRows), x(numCols), res(numRows), Ax(numRows);
      int row_count=1;
      
      for (const_set_it row = testbasis_set.begin(); row != testbasis_set.end(); ++row, ++row_count) {
        rhs(row_count) = truth_model->repr_rhs_F_op((*row));
      }
      timer2.stop();
      std::cout << "    .... done : " << timer2.elapsed() << " seconds " << std::endl;  
      
      T residual;
      int its;
      switch(solution_method){
        case call_cg:
          std::cout << "  Start CG Solve: Maximal iterations = " << maxIterations << std::endl; 
          its = lawa::cg(truth_model->inner_product_matrix, x, rhs, tol, maxIterations);
          Ax = truth_model->inner_product_matrix*x;
          res= Ax-rhs;
          residual = std::sqrt(res*res);
          row_count = 1;
          for (const_set_it row=testbasis_set.begin(); row!=testbasis_set.end(); ++row, ++row_count) {
              u[*row] = x(row_count);
          }
          std::cout << "  CG iterations: " << its << ", residual = " << residual << std::endl;
          break;
        case call_gmres:
          std::cout << "  Start GMRES Solve: Maximal iterations = " << maxIterations << std::endl; 
          its = lawa::gmres(truth_model->inner_product_matrix, x, rhs, tol, maxIterations);
          Ax = truth_model->inner_product_matrix*x;
          res= Ax-rhs;
          residual = std::sqrt(res*res);
          row_count = 1;
          for (const_set_it row=testbasis_set.begin(); row!=testbasis_set.end(); ++row, ++row_count) {
              u[*row] = x(row_count);
          }
          std::cout << "  GMRES iterations: " << its << ", residual = " << residual << std::endl;
          break;
        case call_cgls:
          std::cout << "  Start CGLS Solve: Maximal iterations = " << maxIterations << std::endl; 
          its = lawa::cgls(truth_model->inner_product_matrix, x, rhs, tol, maxIterations);
          Ax = truth_model->inner_product_matrix*x;
          res= Ax-rhs;
          residual = std::sqrt(res*res);
          row_count = 1;
          for (const_set_it row=testbasis_set.begin(); row!=testbasis_set.end(); ++row, ++row_count) {
             u[*row] = x(row_count);
          }
          std::cout << "  CGLS iterations: " << its << ", residual = " << res << std::endl;
          break;
        default: 
          std::cerr << "Method not implemented yet " << std::endl;
          break;
      }
      timer1.stop();
      std::cout << std::endl << "  Done F Representor Solve: "<< timer1.elapsed() << " seconds" << std::endl << std::endl;
      
    }
  }
  
  return u;
}

template <typename T, typename TrialBasis, typename Index, typename Compression, typename TestBasis>
Coefficients<Lexicographical,T,Index>
IndexsetTruthSolver_PG<T, TrialBasis, Index, Compression, TestBasis>::repr_solve_A()
{
  Timer timer1, timer2;
  Coefficients<Lexicographical,T,Index> u;
  
  if(!truth_model->use_inner_product_matrix){ // Assemble LHS Matrix
      
    timer1.start();    
    
    Coefficients<Lexicographical,T,Index> f;  
    // Construct rhs coefficients vector
    typename IndexSet<Index>::const_iterator it;
    typedef typename Coefficients<Lexicographical,T,Index2D>::value_type val_type;
    std::cout << "  Build Rhs... " <<  std::endl;
    timer2.start();
    for(it = basis_set.begin(); it != basis_set.end(); ++it){
      f.insert(val_type((*it), truth_model->repr_rhs_A_op(*it)));
    }
    timer2.stop();
    std::cout << "  .... done : " << timer2.elapsed() << " seconds" << std::endl;

    T res;
    int its;
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
    timer1.stop();
    std::cout << "Done A Representor Solve: "<< timer1.elapsed() << " seconds" << std::endl << std::endl;    
  }
  else{ // Use pre-assembled matrix (assumes assemble_inner_product_matrix has been called before)
    
    std::cout << "  Using inner product matrix!" << std::endl;
    if (basis_set.size() > 0) {
      
      timer1.start();
      
      typedef typename IndexSet<Index>::const_iterator const_set_it;
      typedef typename Coefficients<Lexicographical,T,Index >::const_iterator const_coeff_it;
          
      std::cout << "    Build Dense Vectors ..." << std::endl;
      
      timer2.start();
      int N = basis_set.size();
      DenseVector<Array<T> > rhs(N), x(N), res(N), Ax(N);
      int row_count=1;
      
      for (const_set_it row = basis_set.begin(); row != basis_set.end(); ++row, ++row_count) {
        rhs(row_count) = truth_model->repr_rhs_A_op((*row));
        x(row_count) = 0.;  
      }
      timer2.stop();
      std::cout << "    .... done : " << timer2.elapsed() << " seconds " << std::endl; 
      
      T residual;
      int its;
      switch(solution_method){
        case call_cg:
          std::cout << "  Start CG Solve: Maximal iterations = " << maxIterations << std::endl; 
          its = lawa::cg(truth_model->inner_product_matrix, x, rhs, tol, maxIterations);
          Ax = truth_model->inner_product_matrix*x;
          res= Ax-rhs;
          residual = std::sqrt(res*res);
          row_count = 1;
          for (const_set_it row=basis_set.begin(); row!=basis_set.end(); ++row, ++row_count) {
              u[*row] = x(row_count);
          }
          std::cout << "  CG iterations: " << its << ", residual = " << residual << std::endl;
          break;
        case call_gmres:
          std::cout << "  Start GMRES Solve: Maximal iterations = " << maxIterations << std::endl; 
          its = lawa::gmres(truth_model->inner_product_matrix, x, rhs, tol, maxIterations);
          Ax = truth_model->inner_product_matrix*x;
          res= Ax-rhs;
          residual = std::sqrt(res*res);
          row_count = 1;
          for (const_set_it row=basis_set.begin(); row!=basis_set.end(); ++row, ++row_count) {
              u[*row] = x(row_count);
          }
          std::cout << "  GMRES iterations: " << its << ", residual = " << residual << std::endl;          
          break;
        default: 
          std::cerr << "Method not implemented yet " << std::endl;
          break;
      }
      timer1.stop();
      std::cout << std::endl << "  Done A Representor Solve: "<< timer1.elapsed() << " seconds" << std::endl << std::endl;
      
    }
  }

  
  return u;  
}


template <typename T, typename TrialBasis, typename Index, typename Compression, typename TestBasis>
Coefficients<Lexicographical,T,Index>
IndexsetTruthSolver_PG<T, TrialBasis, Index, Compression, TestBasis>::repr_solve_totalRes(RHS_ResRepr& res_repr_op)
{
  Timer timer1, timer2;
  Coefficients<Lexicographical,T,Index> u;
  
  if(!truth_model->use_inner_product_matrix){ // Assemble LHS Matrix
      
    timer1.start();    
    
    Coefficients<Lexicographical,T,Index> f;  
    // Construct rhs coefficients vector
    typename IndexSet<Index>::const_iterator it;
    typedef typename Coefficients<Lexicographical,T,Index2D>::value_type val_type;
    std::cout << "  Build Rhs... " <<  std::endl;
    timer2.start();
    for(it = basis_set.begin(); it != basis_set.end(); ++it){
      f.insert(val_type((*it), res_repr_op(*it)));
    }
    timer2.stop();
    std::cout << "  .... done : " << timer2.elapsed() << " seconds" << std::endl;

    T res;
    int its;
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
    timer1.stop();
    std::cout << "Done Residual Representor Solve: "<< timer1.elapsed() << " seconds" << std::endl << std::endl;    
  }
  else{ // Use pre-assembled matrix (assumes assemble_inner_product_matrix has been called before)
    
    std::cout << "  Using inner product matrix!" << std::endl;
    if (basis_set.size() > 0) {
      
      timer1.start();
      
      typedef typename IndexSet<Index>::const_iterator const_set_it;
      typedef typename Coefficients<Lexicographical,T,Index >::const_iterator const_coeff_it;
          
      std::cout << "    Build Dense Vectors ..." << std::endl;
      
      timer2.start();
      int N = basis_set.size();
      DenseVector<Array<T> > rhs(N), x(N), res(N), Ax(N);
      int row_count=1;
      
      for (const_set_it row = basis_set.begin(); row != basis_set.end(); ++row, ++row_count) {
        rhs(row_count) = res_repr_op((*row));
        x(row_count) = 0.;  
      }
      timer2.stop();
      std::cout << "    .... done : " << timer2.elapsed() << " seconds " << std::endl; 
      
      T residual;
      int its;
      switch(solution_method){
        case call_cg:
          std::cout << "  Start CG Solve: Maximal iterations = " << maxIterations << std::endl; 
          its = lawa::cg(truth_model->inner_product_matrix, x, rhs, tol, maxIterations);
          Ax = truth_model->inner_product_matrix*x;
          res= Ax-rhs;
          residual = std::sqrt(res*res);
          row_count = 1;
          for (const_set_it row=basis_set.begin(); row!=basis_set.end(); ++row, ++row_count) {
              u[*row] = x(row_count);
          }
          std::cout << "  CG iterations: " << its << ", residual = " << residual << std::endl;
          break;
        case call_gmres:
          std::cout << "  Start GMRES Solve: Maximal iterations = " << maxIterations << std::endl; 
          its = lawa::gmres(truth_model->inner_product_matrix, x, rhs, tol, maxIterations);
          Ax = truth_model->inner_product_matrix*x;
          res= Ax-rhs;
          residual = std::sqrt(res*res);
          row_count = 1;
          for (const_set_it row=basis_set.begin(); row!=basis_set.end(); ++row, ++row_count) {
              u[*row] = x(row_count);
          }
          std::cout << "  GMRES iterations: " << its << ", residual = " << residual << std::endl;          
          break;
        default: 
          std::cerr << "Method not implemented yet " << std::endl;
          break;
      }
      timer1.stop();
      std::cout << std::endl << "  Done Residual Representor Solve: "<< timer1.elapsed() << " seconds" << std::endl << std::endl;
      
    }
  }

  
  return u;
}

} // namespace lawa