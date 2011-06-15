namespace lawa {

template <typename T, typename Basis, typename Index, typename Compression>
IndexsetTruthSolver<T, Basis, Index, Compression>::
IndexsetTruthSolver(IndexSet<Index>& _indexset, Truth& _truth, SolverCall solmethod, T _tol, int maxIts)
  : basis_set(_indexset), truth_model(&_truth), solution_method(solmethod), tol(_tol), maxIterations(maxIts)
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
  
  return u;
}

template <typename T, typename Basis, typename Index, typename Compression>
Coefficients<Lexicographical,T,Index>
IndexsetTruthSolver<T, Basis, Index, Compression>::repr_solve_A()
{
  Coefficients<Lexicographical,T,Index> u, f;
  Timer timer;
  
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
  
  return u;  
}



} // namespace lawa
