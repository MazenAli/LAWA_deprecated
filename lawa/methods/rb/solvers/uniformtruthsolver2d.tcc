namespace lawa {

template <typename T, typename Basis, typename Prec>
UniformTruthSolver2D<T, Basis, Prec>::UniformTruthSolver2D(Basis& _basis)
	: noprec(), prec(noprec), basis(_basis), J_x(_basis.first.j0), J_y(_basis.second.j0), assembler(_basis)
{}


template <typename T, typename Basis, typename Prec>
UniformTruthSolver2D<T, Basis, Prec>::UniformTruthSolver2D(Basis& _basis, int _J_x, int _J_y)
	: noprec(), prec(noprec), basis(_basis), J_x(_J_x), J_y(_J_y), assembler(_basis)
{}

template <typename T, typename Basis, typename Prec>
UniformTruthSolver2D<T, Basis, Prec>::UniformTruthSolver2D(Basis& _basis, Prec& _prec)
	: noprec(), prec(_prec), basis(_basis), assembler(_basis)
{}
        
template <typename T, typename Basis, typename Prec>
UniformTruthSolver2D<T, Basis, Prec>::UniformTruthSolver2D(Basis& _basis, int _J_x, int _J_y, Prec& _prec)
	: noprec(), prec(_prec), basis(_basis), J_x(_J_x), J_y(_J_y), assembler(_basis)
{}


template <typename T, typename Basis, typename Prec>
Coefficients<Lexicographical,T,Index2D>
UniformTruthSolver2D<T, Basis, Prec>::truth_solve()
{
	typedef flens::SparseGeMatrix<flens::CRS<T,flens::CRS_General> > 	SparseMatrixT;
    typedef flens::DenseVector<flens::Array<T> > 						DenseVectorT;
    typedef flens::DiagonalMatrix<T>  									DiagMatrixT;
    
    std::cout << "Assemble System .... Dim = " << basis.dim(J_x, J_y) << std::endl;
    SparseMatrixT A = assembler.assembleStiffnessMatrix(truth_model->lhs_op, J_x, J_y);
    std::cout << ".... A done " << std::endl;
    DenseVectorT  F = assembler.assembleRHS(truth_model->rhs_op, J_x, J_y);
    std::cout << ".... F done " << std::endl;
    DiagMatrixT	  P = assembler.assemblePreconditioner(prec, J_x, J_y);
    std::cout << ".... P done " << std::endl;
    
    std::cout << "Solve System ...." << std::endl;
	DenseVectorT arg(basis.dim(J_x, J_y));
    pcg(P, A, arg, F);
    std::cout << ".... done " << std::endl;

    // VectorToCoefficients
    std::cout << "Convert Solution .... " << std::endl;
    Coefficients<Lexicographical, T, Index2D> u_coeffs;
    denseVectorToCoefficients(arg, u_coeffs);
    std::cout << ".... done " << std::endl;

    return u_coeffs;
}

template <typename T, typename Basis, typename Prec>
void
UniformTruthSolver2D<T, Basis, Prec>::set_Jmax(int _J_x, int _J_y)
{
	J_x = _J_x;
    J_y = _J_y;
}

template <typename T, typename Basis, typename Prec>
void 
UniformTruthSolver2D<T, Basis, Prec>::set_model(Truth& _truth_model){
	truth_model = &_truth_model;
}


template <typename T, typename Basis, typename Prec>
void
UniformTruthSolver2D<T, Basis, Prec>::denseVectorToCoefficients(flens::DenseVector<flens::Array<T> > arg, 
                                                          Coefficients<Lexicographical,T, Index2D>& dest)
{
	dest.clear();
    
	typename Basis::FirstBasisType b1 = basis.first;
    typename Basis::SecondBasisType b2 = basis.second;
    typedef typename Coefficients<Lexicographical,T,Index2D>::value_type val_type;

    UniformIndex2D<Basis> I(basis, J_x, J_y);
   
     /* We have to convert each basis function into an index 
     * and insert that and its corresponding value in arg into dest
     */
    
    // ScalingFct x ScalingFct
    Range<int> Rux = b1.mra.rangeI(b1.j0);
    Range<int> Ruy = b2.mra.rangeI(b2.j0);
    for (int kux = Rux.firstIndex(); kux <= Rux.lastIndex(); kux++) {
        for (int kuy = Ruy.firstIndex(); kuy <= Ruy.lastIndex(); kuy++) {
        	Index1D ind_x(b1.j0, kux, XBSpline);
            Index1D ind_y(b2.j0, kuy, XBSpline);
            Index2D ind_2d(ind_x, ind_y);
            dest.insert(val_type(ind_2d, arg(I(XBSpline, b1.j0, kux, XBSpline, b2.j0, kuy))));
        }
    }
    
    // ScalingFct x Wavelet
    Rux = b1.mra.rangeI(b1.j0);
    for (int juy = b2.j0; juy < basis.J2_max(J_x, J_y, b1.j0-1); juy++) {
        Ruy = b2.rangeJ(juy);
        for(int kux = Rux.firstIndex(); kux <= Rux.lastIndex(); ++kux){
    		for(int kuy = Ruy.firstIndex(); kuy <= Ruy.lastIndex(); ++kuy){       
                Index1D ind_x(b1.j0, kux, XBSpline);
                Index1D ind_y(  juy, kuy, XWavelet);
                Index2D ind_2d(ind_x, ind_y);
                dest.insert(val_type(ind_2d, arg(I(XBSpline, b1.j0, kux, XWavelet, juy, kuy))));            	
            }
        }
    }
    
    // Wavelet x ScalingFct
    Ruy = b2.mra.rangeI(b2.j0);
    for (int jux = b1.j0; jux < basis.J1_max(J_x, J_y, b2.j0-1); jux++) {
    	Rux = b1.rangeJ(jux);
        for(int kux = Rux.firstIndex(); kux <= Rux.lastIndex(); ++kux){
    		for(int kuy = Ruy.firstIndex(); kuy <= Ruy.lastIndex(); ++kuy){       
                Index1D ind_x(	jux, kux, XWavelet);
                Index1D ind_y(b2.j0, kuy, XBSpline);
                Index2D ind_2d(ind_x, ind_y);
                dest.insert(val_type(ind_2d, arg(I(XWavelet, jux, kux, XBSpline, b2.j0, kuy))));            	
            }
        }        
    }
    
    // Wavelet x Wavelet
    for(int jux = b1.j0; jux < basis.J1_max(J_x, J_y, b2.j0); ++jux){
    	Rux = b1.rangeJ(jux); 
       	for(int juy = b2.j0; juy < basis.J2_max(J_x, J_y, jux); ++juy){    
        	Ruy = b2.rangeJ(juy); 
         	for(int kux = Rux.firstIndex(); kux <= Rux.lastIndex(); ++kux){
           		for(int kuy = Ruy.firstIndex(); kuy <= Ruy.lastIndex(); ++kuy){
                    Index1D ind_x(jux, kux, XWavelet);
                    Index1D ind_y(juy, kuy, XWavelet);
                    Index2D ind_2d(ind_x, ind_y);
                    dest.insert(val_type(ind_2d, arg(I(XWavelet, jux, kux, XWavelet, juy, kuy)))); 
           		}
        	}
    	}            
    }    
    
}

} // namespace lawa
