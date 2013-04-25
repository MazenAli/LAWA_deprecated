#include <fstream>

namespace lawa {

template <typename T, typename Basis, typename Prec, typename Index, typename Compression>
IndexsetTruthSolver<T, Basis, Prec, Index, Compression>::
IndexsetTruthSolver(IndexSet<Index>& _indexset, Truth& _truth, SolverCall solmethod, T _tol, int maxIts)
  : basis_set(_indexset), truth_model(&_truth), solution_method(solmethod), tol(_tol), maxIterations(maxIts)
{}

template <typename T, typename Basis, typename Prec, typename Index, typename Compression>
void 
IndexsetTruthSolver<T, Basis, Prec, Index, Compression>
::set_model(AdaptiveRBTruth2D<T, Basis, Prec, IndexsetTruthSolver<T, Basis, Prec, Index, Compression>, Compression>& _truth_model){
    truth_model = &_truth_model;
}

template <typename T, typename Basis, typename Prec, typename Index, typename Compression>
Coefficients<Lexicographical,T,Index>
IndexsetTruthSolver<T, Basis, Prec, Index, Compression>::truth_solve()
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
  for(it = basis_set.begin(); it != basis_set.end(); ++it){
    f.insert(val_type((*it), truth_model->rhs_op(*it)));
  }
  timer2.stop();
  std::cout << "  ..... done : " << timer2.elapsed() << " seconds" << std::endl;
  
  T res;
  int its;
	int assemble_matrix = 1;
  switch(solution_method){
    case call_cg:
		{
			T timeMatrixVector=0;
    	std::cout << "  Start CG Solve: Maximal iterations = " << maxIterations << std::endl; 
      its = CG_Solve(basis_set, truth_model->lhs_op, u, f, res, tol, maxIterations, timeMatrixVector, assemble_matrix);
      std::cout << "  CG iterations: " << its << ", residual = " << res << std::endl;
      break;
		}
    case call_gmres:
    std::cout << "  Start GMRES Solve: Maximal iterations = " << maxIterations << std::endl; 
      its = GMRES_Solve(basis_set, truth_model->lhs_op, u, f, res, tol, maxIterations, assemble_matrix);
      std::cout << "  GMRES iterations: " << its << ", residual = " << res << std::endl;
      break;
    default: 
      std::cerr << "Method not implemented yet " << std::endl;
      break;
  }
  timer1.stop();
  std::cout << "Done Truth Solve: "<< timer1.elapsed() << " seconds" << std::endl << std::endl;
  
  return u;
}

template <typename T, typename Basis, typename Prec, typename Index, typename Compression>
Coefficients<Lexicographical,T,Index>
IndexsetTruthSolver<T, Basis, Prec, Index, Compression>::repr_solve_F()
{
  Coefficients<Lexicographical,T,Index> u, f;
  Timer timer1, timer2;
  
  timer1.start();
  
  if(!truth_model->use_inner_product_matrix){ // Assemble LHS Matrix
  	int assemble_matrix = 1;
		T timeMatrixVector=0;
	
    // Construct rhs coefficients vector
    typename IndexSet<Index>::const_iterator it;
    typedef typename Coefficients<Lexicographical,T,Index2D>::value_type val_type;
    std::cout << "  Build Rhs... " <<  std::endl;
    timer2.start();
    
    for(it = basis_set.begin(); it != basis_set.end(); ++it){
      f.insert(val_type((*it), truth_model->repr_rhs_F_op(*it)));
    }
    timer2.stop();
    std::cout << "  .... done : " << timer2.elapsed() << " seconds" << std::endl;
  
    T res;
    int its;

    std::cout << "  Start CG Solve: Maximal iterations = " << maxIterations << std::endl; 
    its = CG_Solve(basis_set, truth_model->repr_lhs_op, u, f, res, tol, maxIterations, timeMatrixVector, assemble_matrix);
    std::cout << "  CG iterations: " << its << ", residual = " << res << std::endl;

    timer1.stop();
    std::cout << "Done F Representor Solve: "<< timer1.elapsed() << " seconds" << std::endl << std::endl;
  }
  else{ // Use pre-assembled matrix (assumes assemble_inner_product_matrix has been called before)
    
    assert(truth_model->assembled_inner_product_matrix);
    assert(truth_model->assembled_prec_vec);
    
    DiagonalMatrix<T> P(truth_model->prec_vec);
    
    std::cout << "    Using inner product matrix!" << std::endl;
    if (basis_set.size() > 0) {
      
      typedef typename IndexSet<Index>::const_iterator const_set_it;
      typedef typename Coefficients<Lexicographical,T,Index >::const_iterator const_coeff_it;
          
      std::cout << "    Build Dense Vectors ..." << std::endl;
      
      timer2.start();
      int N = basis_set.size();
      DenseVector<Array<T> > rhs(N), x(N), res(N), Ax(N);
      int row_count=1;

      for (const_set_it row = basis_set.begin(); row != basis_set.end(); ++row, ++row_count) {
         rhs(row_count) = truth_model->repr_rhs_F_op((*row)) / truth_model->get_prec((*row));

        x(row_count) = 0.;  
      }
      timer2.stop();
      std::cout << "    .... done : " << timer2.elapsed() << " seconds " << std::endl;  
      
      T residual;
      int its;

      std::cout << "  Start PCG Solve: Maximal iterations = " << maxIterations << std::endl; 
      its = lawa::pcg(P, truth_model->inner_product_matrix, x, rhs, tol, maxIterations);
      Ax = truth_model->inner_product_matrix*x;
      res= Ax-rhs;
      residual = std::sqrt(res*res);
      row_count = 1;
      for (const_set_it row=basis_set.begin(); row!=basis_set.end(); ++row, ++row_count) {
          u[*row] = x(row_count) / truth_model->get_prec((*row)) ;
      }
      std::cout << "  PCG iterations: " << its << ", residual = " << residual << std::endl;

      timer1.stop();
      std::cout << std::endl << "  Done F Representor Solve: "<< timer1.elapsed() << " seconds" << std::endl << std::endl;
      
    }
  }
  
  return u;
}

template <typename T, typename Basis, typename Prec, typename Index, typename Compression>
Coefficients<Lexicographical,T,Index>
IndexsetTruthSolver<T, Basis, Prec, Index, Compression>::repr_solve_A()
{
  Timer timer1, timer2;
  Coefficients<Lexicographical,T,Index> u;
  
  if(!truth_model->use_inner_product_matrix){ // Assemble LHS Matrix
		
		int assemble_matrix = 1;
		T timeMatrixVector=0;
		
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

    std::cout << "  Start CG Solve: Maximal iterations = " << maxIterations << std::endl; 
    its = CG_Solve(basis_set, truth_model->repr_lhs_op, u, f, res, tol, maxIterations, timeMatrixVector, assemble_matrix);
    std::cout << "  CG iterations: " << its << ", residual = " << res << std::endl;

    timer1.stop();
    std::cout << "Done A Representor Solve: "<< timer1.elapsed() << " seconds" << std::endl << std::endl;    
  }
  else{ // Use pre-assembled matrix (assumes assemble_inner_product_matrix has been called before)
    
    assert(truth_model->assembled_inner_product_matrix);
    assert(truth_model->assembled_prec_vec);
      
    DiagonalMatrix<T> P(truth_model->prec_vec);  
    
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
        rhs(row_count) = truth_model->repr_rhs_A_op((*row)) / truth_model->get_prec((*row));
        x(row_count) = 0.;  
      }
      timer2.stop();
      std::cout << "    .... done : " << timer2.elapsed() << " seconds " << std::endl; 
      
      T residual;
      int its;

      std::cout << "  Start PCG Solve: Maximal iterations = " << maxIterations << std::endl; 
      its = lawa::pcg(P, truth_model->inner_product_matrix, x, rhs, tol, maxIterations);
      Ax = truth_model->inner_product_matrix*x;
      res= Ax-rhs;
      residual = std::sqrt(res*res);
      row_count = 1;
      for (const_set_it row=basis_set.begin(); row!=basis_set.end(); ++row, ++row_count) {
          u[*row] = x(row_count) / truth_model->get_prec((*row));
      }
      std::cout << "  PCG iterations: " << its << ", residual = " << residual << std::endl;
      
      timer1.stop();
      std::cout << std::endl << "  Done A Representor Solve: "<< timer1.elapsed() << " seconds" << std::endl << std::endl;
      
    }
  }

  
  return u;  
}

template <typename T, typename Basis, typename Prec, typename Index, typename Compression>
Coefficients<Lexicographical,T,Index>
IndexsetTruthSolver<T, Basis, Prec, Index, Compression>::repr_solve_output()
{
  Coefficients<Lexicographical,T,Index> u, f;
  Timer timer1, timer2;

  timer1.start();

  if(!truth_model->use_inner_product_matrix){ // Assemble LHS Matrix
   
		int assemble_matrix = 1;
		T timeMatrixVector=0;
		
        // Construct rhs coefficients vector
    typename IndexSet<Index>::const_iterator it;
    typedef typename Coefficients<Lexicographical,T,Index2D>::value_type val_type;
    std::cout << "  Build Rhs... " <<  std::endl;
    timer2.start();
    for(it = basis_set.begin(); it != basis_set.end(); ++it){
      f.insert(val_type((*it), truth_model->repr_rhs_output_op(*it)));
    }
    timer2.stop();
    std::cout << "  .... done : " << timer2.elapsed() << " seconds" << std::endl;

    T res;
    int its;

    std::cout << "  Start CG Solve: Maximal iterations = " << maxIterations << std::endl; 
    its = CG_Solve(basis_set, truth_model->repr_lhs_op, u, f, res, tol, maxIterations, timeMatrixVector, assemble_matrix);
    std::cout << "  CG iterations: " << its << ", residual = " << res << std::endl;

    timer1.stop();
    std::cout << "Done Output Representor Solve: "<< timer1.elapsed() << " seconds" << std::endl << std::endl;
  }
  else{ // Use pre-assembled matrix (assumes assemble_inner_product_matrix has been called before)

    assert(truth_model->assembled_inner_product_matrix);
    assert(truth_model->assembled_prec_vec);

		DiagonalMatrix<T> P(truth_model->prec_vec);  

    std::cout << "    Using inner product matrix!" << std::endl;
    if (basis_set.size() > 0) {

      typedef typename IndexSet<Index>::const_iterator const_set_it;
      typedef typename Coefficients<Lexicographical,T,Index >::const_iterator const_coeff_it;

      std::cout << "    Build Dense Vectors ..." << std::endl;
      
      timer2.start();
      int N = basis_set.size();
      DenseVector<Array<T> > rhs(N), x(N), res(N), Ax(N);
      int row_count=1;

      for (const_set_it row = basis_set.begin(); row != basis_set.end(); ++row, ++row_count) {
        rhs(row_count) = truth_model->repr_rhs_output_op((*row))/ truth_model->get_prec((*row));
        x(row_count) = 0.;
      }
      timer2.stop();
      std::cout << "    .... done : " << timer2.elapsed() << " seconds " << std::endl;

      T residual;
      int its;
      std::cout << "  Start PCG Solve: Maximal iterations = " << maxIterations << std::endl; 
      its = lawa::pcg(P, truth_model->inner_product_matrix, x, rhs, tol, maxIterations);
      Ax = truth_model->inner_product_matrix*x;
      res= Ax-rhs;
      residual = std::sqrt(res*res);
      row_count = 1;
      for (const_set_it row=basis_set.begin(); row!=basis_set.end(); ++row, ++row_count) {
          u[*row] = x(row_count) / truth_model->get_prec((*row));
      }
      std::cout << "  PCG iterations: " << its << ", residual = " << residual << std::endl;

      timer1.stop();
      std::cout << std::endl << "  Done Output Representor Solve: "<< timer1.elapsed() << " seconds" << std::endl << std::endl;

    }
  }

  return u;
}


template <typename T, typename Basis, typename Prec, typename Index, typename Compression>
Coefficients<Lexicographical,T,Index>
IndexsetTruthSolver<T, Basis, Prec, Index, Compression>::repr_solve_totalRes(RHS_ResRepr& res_repr_op)
{
  Timer timer1, timer2;
  Coefficients<Lexicographical,T,Index> u;
  
  if(!truth_model->use_inner_product_matrix){ // Assemble LHS Matrix
      
    timer1.start();    

		int assemble_matrix = 1;
		T timeMatrixVector = 0;
    
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
    std::cout << "  Start CG Solve: Maximal iterations = " << maxIterations << std::endl; 
    its = CG_Solve(basis_set, truth_model->repr_lhs_op, u, f, res, tol, maxIterations, timeMatrixVector, assemble_matrix);
    std::cout << "  CG iterations: " << its << ", residual = " << res << std::endl;

    timer1.stop();
    std::cout << "Done Residual Representor Solve: "<< timer1.elapsed() << " seconds" << std::endl << std::endl;    
  }
  else{ // Use pre-assembled matrix (assumes assemble_inner_product_matrix has been called before)
    
    assert(truth_model->assembled_inner_product_matrix);
    assert(truth_model->assembled_prec_vec);
      
    DiagonalMatrix<T> P(truth_model->prec_vec);  
    
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
        rhs(row_count) = res_repr_op((*row)) / truth_model->get_prec((*row));
        x(row_count) = 0.;  
      }
      timer2.stop();
      std::cout << "    .... done : " << timer2.elapsed() << " seconds " << std::endl; 
      
      T residual;
      int its;

      std::cout << "  Start PCG Solve: Maximal iterations = " << maxIterations << std::endl; 
      its = lawa::pcg(P, truth_model->inner_product_matrix, x, rhs, tol, maxIterations);
      Ax = truth_model->inner_product_matrix*x;
      res= Ax-rhs;
      residual = std::sqrt(res*res);
      row_count = 1;
      for (const_set_it row=basis_set.begin(); row!=basis_set.end(); ++row, ++row_count) {
          u[*row] = x(row_count) / truth_model->get_prec((*row));
      }
      std::cout << "  PCG iterations: " << its << ", residual = " << residual << std::endl;

      timer1.stop();
      std::cout << std::endl << "  Done Residual Representor Solve: "<< timer1.elapsed() << " seconds" << std::endl << std::endl;
      
    }
  }

  return u;
}

} // namespace lawa