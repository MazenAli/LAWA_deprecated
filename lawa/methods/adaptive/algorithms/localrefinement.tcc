namespace lawa {

template <typename PrimalBasis>
LocalRefinement<PrimalBasis>::LocalRefinement(const PrimalBasis &_basis, bool withDirichletBC)
 : basis(_basis),
   numLeftBoundaryScalings(basis.mra.cardIL(basis.j0)),
   numRightBoundaryScalings(basis.mra.cardIR(basis.j0)),
   numLeftBoundaryWavelets(basis.cardJL(basis.j0)),
   numRightBoundaryWavelets(basis.cardJR(basis.j0)),
   inner_scaling_offset(0), inner_wavelet_offset(0),
   leftM0(numLeftBoundaryScalings), rightM0(numRightBoundaryScalings)
{
    const int d =basis.d;
    const int d_=basis.d_;
    const int j = ((d==2) && ((d_==2)||(d_==4))) ? basis.mra_.min_j0+1 : basis.mra_.min_j0;

    GeMatrix<FullStorage<T,ColMajor> > Mj1, Mj1_;
    initial_stable_completion(basis.mra,basis.mra_,Mj1,Mj1_);

    if (Mj1.numCols()!=basis.cardJ(j) && Mj1.numRows()!=basis.mra.cardI(j+1)) {
        std::cerr << "LocalRefinement: Refinement matrix does not have correct dimension."
                  << std::endl;
    }

    innerM0.engine().resize(basis.mra.M0.rightband.length());
    innerM0(_(1,basis.mra.M0.rightband.length())) = basis.mra.M0.rightband;

    leftM0.engine().resize(basis.mra.M0.left.length());
    for (int i=0; i<basis.mra.M0.left.length(); ++i) {
        int pos = basis.mra.M0.left.firstIndex()+i;
        DenseVector<Array<T> > left(basis.mra.M0.left(pos).length());
        left(_(1,basis.mra.M0.left(pos).length())) = basis.mra.M0.left(pos);
        leftM0(i+1) = left;
    }

    rightM0.engine().resize(basis.mra.M0.left.length());
    for (int i=0; i<basis.mra.M0.right.length(); ++i) {
        int pos = basis.mra.M0.right.lastIndex()-i;
        DenseVector<Array<T> > right(basis.mra.M0.right(pos).length());
        right(_(1,basis.mra.M0.right(pos).length())) = basis.mra.M0.right(pos);
        rightM0(i+1) = right;
    }


    if (basis.Cons==Dijkema) {
        if (d==2 && d_==2) {
            inner_wavelet_offset = 4;
            if (withDirichletBC) {  inner_scaling_offset = 0;  }
            else                 {  inner_scaling_offset = 2;  }
        }
        else if (d==3 && d_==3) {
            inner_wavelet_offset = 7;
            if (withDirichletBC) {  inner_scaling_offset = 1;   }
            else                 {  inner_scaling_offset = 3;   }
        }
    }

    innerM1.engine().resize(basis.M1.leftband.length());
    innerM1(_(1,basis.M1.leftband.length())) = basis.M1.leftband;

    leftM1.engine().resize(numLeftBoundaryWavelets);
    for (int i=0; i<numLeftBoundaryWavelets; ++i) {
        int pos = basis.M1.left.firstIndex()+i;
        DenseVector<Array<T> > left(basis.M1.left(pos).length());
        left(_(1,basis.M1.left(pos).length())) = basis.M1.left(pos);
        leftM1(i+1) = left;
    }

    rightM1.engine().resize(numRightBoundaryWavelets);
    for (int i=0; i<numLeftBoundaryWavelets; ++i) {
        int pos = basis.M1.right.lastIndex()-i;
        DenseVector<Array<T> > right(basis.M1.right(pos).length());
        right(_(1,basis.M1.right(pos).length())) = basis.M1.right(pos);
        rightM1(i+1) = right;
    }
}

template <typename PrimalBasis>
void
LocalRefinement<PrimalBasis>::reconstruct(const Coefficients<Lexicographical,T,Index1D> &u_multi_j,
                                          Coefficients<Lexicographical,T,Index1D> &u_loc_single_jP1) const
{
    for (const_coeff1d_it it=u_multi_j.begin(); it!=u_multi_j.end(); ++it) {
        this->reconstruct((*it).first.j, (*it).first.k, (*it).first.xtype,
                          (*it).second, u_loc_single_jP1);
    }
}

template <typename PrimalBasis>
void
LocalRefinement<PrimalBasis>::reconstruct(const Coefficients<Lexicographical,T,Index1D> &u_multi_j,
                                          int j,
                                          Coefficients<Lexicographical,T,Index1D> &u_loc_single_jP1) const
{
    for (const_coeff1d_it it=u_multi_j.begin(); it!=u_multi_j.end(); ++it) {
        if ((*it).first.j==j) {
            this->reconstruct((*it).first.j, (*it).first.k, (*it).first.xtype,
                              (*it).second, u_loc_single_jP1);
        }
    }
}

template <typename PrimalBasis>
void
LocalRefinement<PrimalBasis>::reconstruct(const CoefficientsByLevel<T> &u_scaling,
                                          const CoefficientsByLevel<T> &u_wavelet, int j,
                                          CoefficientsByLevel<T> &u_loc_single_jP1) const
{
    for (typename CoefficientsByLevel<T>::const_it it=u_scaling.map.begin(); it!=u_scaling.map.end(); ++it) {
        this->reconstruct(j, (*it).first, XBSpline, (*it).second, u_loc_single_jP1);
    }
    for (typename CoefficientsByLevel<T>::const_it it=u_wavelet.map.begin(); it!=u_wavelet.map.end(); ++it) {
        this->reconstruct(j, (*it).first, XWavelet, (*it).second, u_loc_single_jP1);
    }
}

template <typename PrimalBasis>
void
LocalRefinement<PrimalBasis>::reconstruct(const short &j, const long &k, const XType &xtype, T coeff,
                                          Coefficients<Lexicographical,T,Index1D> &u_loc_single) const
{
    if (xtype==XWavelet) {
        if (k<=numLeftBoundaryWavelets) {
            long k_first_scaling=basis.mra.rangeI(j+1).firstIndex()-1;
            for (int i=leftM1(k).firstIndex(); i<=leftM1(k).lastIndex(); ++i) {
                u_loc_single[Index1D(j+1,k_first_scaling + i,XBSpline)] += coeff * leftM1(k).operator()(i);
            }
        }
        else if (k>basis.cardJ(j)-numRightBoundaryWavelets) {
            int type = basis.cardJ(j) - k + 1;
            long k_first_scaling= basis.mra.rangeI(j+1).lastIndex() - rightM1(type).length();
            for (int i=rightM1(type).firstIndex(); i<=rightM1(type).lastIndex(); ++i) {
                u_loc_single[Index1D(j+1,k_first_scaling+i,XBSpline)] += coeff * rightM1(type).operator()(i);
            }
        }
        else {
            long k_first_scaling= (k - numLeftBoundaryWavelets - 1)*2 + inner_wavelet_offset - 1;
            T factor = ((basis.d % 2 != 0) && (k>basis.cardJ(j)/2.)) ? -1. : 1.;
            for (int i=innerM1.firstIndex(); i<=innerM1.lastIndex(); ++i) {
                u_loc_single[Index1D(j+1,k_first_scaling+i,XBSpline)] += coeff * factor * innerM1(i);
            }
        }
    }
    else {
        if (k<=basis.mra.rangeIL(j).lastIndex()) {
            int offset = basis.mra.rangeI(basis.j0).firstIndex()-1;
            long k_first_scaling= basis.mra.rangeI(j+1).firstIndex() - 1;
            for (int i=leftM0(k-offset).firstIndex(); i<=leftM0(k-offset).lastIndex(); ++i) {
                u_loc_single[Index1D(j+1,k_first_scaling+i,XBSpline)] += coeff * leftM0(k-offset).operator()(i);
            }
        }
        else if (k>=basis.mra.rangeIR(j).firstIndex()) {
            int type = basis.mra.rangeI(j).lastIndex() - k + 1;
            long k_first_scaling= basis.mra.rangeI(j+1).lastIndex() - rightM0(type).length();
            for (int i=rightM0(type).firstIndex(); i<=rightM0(type).lastIndex(); ++i) {
                u_loc_single[Index1D(j+1,k_first_scaling+i,XBSpline)] += coeff * rightM0(type).operator()(i);
            }
        }
        else {
            long k_first_scaling= (k - numLeftBoundaryScalings - 1)*2 + inner_scaling_offset - 1;
            for (int i=innerM0.firstIndex(); i<=innerM0.lastIndex(); ++i) {
                u_loc_single[Index1D(j+1,k_first_scaling+i,XBSpline)] += coeff * innerM0(i);
            }
        }
    }
}

template <typename PrimalBasis>
void
LocalRefinement<PrimalBasis>::reconstruct(const short &j, const long &k, const XType &xtype, T coeff,
                                          CoefficientsByLevel<T> &u_loc_single) const
{
    if (xtype==XWavelet) {
        if (k<=numLeftBoundaryWavelets) {
            long k_first_scaling=basis.mra.rangeI(j+1).firstIndex()-1;
            for (int i=leftM1(k).firstIndex(); i<=leftM1(k).lastIndex(); ++i) {
                u_loc_single.map[k_first_scaling+i] += coeff * leftM1(k).operator()(i);
            }
        }
        else if (k>basis.cardJ(j)-numRightBoundaryWavelets) {
            int type = basis.cardJ(j) - k + 1;
            long k_first_scaling= basis.mra.rangeI(j+1).lastIndex() - rightM1(type).length();
            for (int i=rightM1(type).firstIndex(); i<=rightM1(type).lastIndex(); ++i) {
                u_loc_single.map[k_first_scaling+i] += coeff * rightM1(type).operator()(i);
            }
        }
        else {
            long k_first_scaling= (k - numLeftBoundaryWavelets - 1)*2 + inner_wavelet_offset - 1;
            T factor = ((basis.d % 2 != 0) && (k>basis.cardJ(j)/2.)) ? -1. : 1.;
            for (int i=innerM1.firstIndex(); i<=innerM1.lastIndex(); ++i) {
                u_loc_single.map[k_first_scaling+i] += coeff * factor * innerM1(i);
            }
        }
    }
    else {
        if (k<=basis.mra.rangeIL(j).lastIndex()) {
            int offset = basis.mra.rangeI(basis.j0).firstIndex()-1;
            long k_first_scaling= basis.mra.rangeI(j+1).firstIndex() - 1;
            for (int i=leftM0(k-offset).firstIndex(); i<=leftM0(k-offset).lastIndex(); ++i) {
                u_loc_single.map[k_first_scaling+i] += coeff * leftM0(k-offset).operator()(i);
            }
        }
        else if (k>=basis.mra.rangeIR(j).firstIndex()) {
            int type = basis.mra.rangeI(j).lastIndex() - k + 1;
            long k_first_scaling= basis.mra.rangeI(j+1).lastIndex() - rightM0(type).length();
            for (int i=rightM0(type).firstIndex(); i<=rightM0(type).lastIndex(); ++i) {
                u_loc_single.map[k_first_scaling+i] += coeff * rightM0(type).operator()(i);
            }
        }
        else {
            long k_first_scaling= (k - numLeftBoundaryScalings - 1)*2 + inner_scaling_offset - 1;
            for (int i=innerM0.firstIndex(); i<=innerM0.lastIndex(); ++i) {
                u_loc_single.map[k_first_scaling+i] += coeff * innerM0(i);
            }
        }
    }
}



template <typename PrimalBasis>
void
LocalRefinement<PrimalBasis>::decompose_(Coefficients<Lexicographical,T,Index1D> &u_loc_single,
                                         const IndexSet<Index1D> &Lambda,
                                         Coefficients<Lexicographical,T,Index1D> &u_loc_single_jM1,
                                         Coefficients<Lexicographical,T,Index1D> &u_multi) const
{
    for (const_set1d_it row=Lambda.begin(); row!=Lambda.end(); ++row) {
        T coeff = this->decompose_(u_loc_single, (*row).j, (*row).k, (*row).xtype);
        if ((*row).xtype==XBSpline) u_loc_single_jM1[*row] += coeff;
        else                        u_multi[*row] += coeff;
    }
}

template <typename PrimalBasis>
void
LocalRefinement<PrimalBasis>::decompose_(CoefficientsByLevel<T>  &u_loc_single, int j,
                                         CoefficientsByLevel<T>  &u_scaling,
                                         CoefficientsByLevel<T>  &u_wavelet) const
{
    if (u_loc_single.map.size()==0) return;
    for (typename CoefficientsByLevel<T>::iter it=u_scaling.map.begin(); it!=u_scaling.map.end(); ++it) {
        T coeff = this->decompose_(u_loc_single, j-1, (*it).first, XBSpline);
        (*it).second += coeff;
    }
    for (typename CoefficientsByLevel<T>::iter it=u_wavelet.map.begin(); it!=u_wavelet.map.end(); ++it) {
        T coeff = this->decompose_(u_loc_single, j-1, (*it).first, XWavelet);
        (*it).second += coeff;
    }
}

template <typename PrimalBasis>
typename PrimalBasis::T
LocalRefinement<PrimalBasis>::decompose_(Coefficients<Lexicographical,T,Index1D> &u_loc_single,
                                         const short &j, const long &k, const XType &xtype) const
{
    T val = 0.;
    if (xtype==XWavelet) {
        if (k<=numLeftBoundaryWavelets) {
            long k_first_scaling=basis.mra.rangeI(j+1).firstIndex()-1;
            for (int i=leftM1(k).firstIndex(); i<=leftM1(k).lastIndex(); ++i) {
                val += u_loc_single[Index1D(j+1,k_first_scaling + i,XBSpline)] * leftM1(k).operator()(i);
            }
        }
        else if (k>basis.cardJ(j)-numRightBoundaryWavelets) {
            int type = basis.cardJ(j) - k + 1;
            long k_first_scaling= basis.mra.rangeI(j+1).lastIndex() - rightM1(type).length();
            for (int i=rightM1(type).firstIndex(); i<=rightM1(type).lastIndex(); ++i) {
                val += u_loc_single[Index1D(j+1,k_first_scaling+i,XBSpline)] * rightM1(type).operator()(i);
            }
        }
        else {
            long k_first_scaling= (k - numLeftBoundaryWavelets - 1)*2 + inner_wavelet_offset - 1;
            T factor = ((basis.d % 2 != 0) && (k>basis.cardJ(j)/2.)) ? -1. : 1.;
            for (int i=innerM1.firstIndex(); i<=innerM1.lastIndex(); ++i) {
                val += u_loc_single[Index1D(j+1,k_first_scaling+i,XBSpline)] * factor * innerM1(i);
            }
        }
    }
    else {
        if (k<=basis.mra.rangeIL(j).lastIndex()) {
            int offset = basis.mra.rangeI(basis.j0).firstIndex()-1;
            long k_first_scaling= basis.mra.rangeI(j+1).firstIndex() - 1;
            for (int i=leftM0(k-offset).firstIndex(); i<=leftM0(k-offset).lastIndex(); ++i) {
                val += u_loc_single[Index1D(j+1,k_first_scaling+i,XBSpline)] * leftM0(k-offset).operator()(i);
            }
        }
        else if (k>=basis.mra.rangeIR(j).firstIndex()) {
            int type = basis.mra.rangeI(j).lastIndex() - k + 1;
            long k_first_scaling= basis.mra.rangeI(j+1).lastIndex() - rightM0(type).length();
            for (int i=rightM0(type).firstIndex(); i<=rightM0(type).lastIndex(); ++i) {
                val += u_loc_single[Index1D(j+1,k_first_scaling+i,XBSpline)] * rightM0(type).operator()(i);
            }
        }
        else {
            long k_first_scaling= (k - numLeftBoundaryScalings - 1)*2 + inner_scaling_offset - 1;
            for (int i=innerM0.firstIndex(); i<=innerM0.lastIndex(); ++i) {
                val += u_loc_single[Index1D(j+1,k_first_scaling+i,XBSpline)] * innerM0(i);
            }
        }
    }
    return val;
}

template <typename PrimalBasis>
typename PrimalBasis::T
LocalRefinement<PrimalBasis>::decompose_(const CoefficientsByLevel<T> &u_loc_single,
                                         const short &j, const long &k, const XType &xtype) const
{
    const_coeffbylevel_it u_loc_single_end = u_loc_single.map.end();
    const_coeffbylevel_it u_loc_single_ptr;
    long double val = 0.L;
    if (xtype==XWavelet) {
        if (k<=numLeftBoundaryWavelets) {
            long k_first_scaling=basis.mra.rangeI(j+1).firstIndex()-1;
            for (int i=leftM1(k).firstIndex(); i<=leftM1(k).lastIndex(); ++i) {
                u_loc_single_ptr=u_loc_single.map.find(k_first_scaling+i);
                if (u_loc_single_ptr!=u_loc_single_end) {
                    val += (*u_loc_single_ptr).second * leftM1(k).operator()(i);
                }
                //val += (long double)(u_loc_single.map[k_first_scaling+i] * leftM1(k).operator()(i));
            }
        }
        else if (k>basis.cardJ(j)-numRightBoundaryWavelets) {
            int type = basis.cardJ(j) - k + 1;
            long k_first_scaling= basis.mra.rangeI(j+1).lastIndex() - rightM1(type).length();
            for (int i=rightM1(type).firstIndex(); i<=rightM1(type).lastIndex(); ++i) {
                u_loc_single_ptr=u_loc_single.map.find(k_first_scaling+i);
                if (u_loc_single_ptr!=u_loc_single_end) {
                    val += (*u_loc_single_ptr).second * rightM1(type).operator()(i);
                }
                //val += (long double)(u_loc_single.map[k_first_scaling+i] * rightM1(type).operator()(i));
            }
        }
        else {
            long k_first_scaling= (k - numLeftBoundaryWavelets - 1)*2 + inner_wavelet_offset - 1;
            T factor = ((basis.d % 2 != 0) && (k>basis.cardJ(j)/2.)) ? -1. : 1.;
            for (int i=innerM1.firstIndex(); i<=innerM1.lastIndex(); ++i) {
                u_loc_single_ptr=u_loc_single.map.find(k_first_scaling+i);
                if (u_loc_single_ptr!=u_loc_single_end) {
                    val += (*u_loc_single_ptr).second * factor * innerM1(i);
                }
                //val += (long double)(u_loc_single.map[k_first_scaling+i] * factor * innerM1(i));
            }
        }
    }
    else {
        if (k<=basis.mra.rangeIL(j).lastIndex()) {
            int offset = basis.mra.rangeI(basis.j0).firstIndex()-1;
            long k_first_scaling= basis.mra.rangeI(j+1).firstIndex() - 1;
            for (int i=leftM0(k-offset).firstIndex(); i<=leftM0(k-offset).lastIndex(); ++i) {
                u_loc_single_ptr=u_loc_single.map.find(k_first_scaling+i);
                if (u_loc_single_ptr!=u_loc_single_end) {
                    val += (*u_loc_single_ptr).second * leftM0(k-offset).operator()(i);
                }
                //val += (long double)(u_loc_single.map[k_first_scaling+i] * leftM0(k-offset).operator()(i));
            }
        }
        else if (k>=basis.mra.rangeIR(j).firstIndex()) {
            int type = basis.mra.rangeI(j).lastIndex() - k + 1;
            long k_first_scaling= basis.mra.rangeI(j+1).lastIndex() - rightM0(type).length();
            for (int i=rightM0(type).firstIndex(); i<=rightM0(type).lastIndex(); ++i) {
                u_loc_single_ptr=u_loc_single.map.find(k_first_scaling+i);
                if (u_loc_single_ptr!=u_loc_single_end) {
                    val += (*u_loc_single_ptr).second * rightM0(type).operator()(i);
                }
                //val += (long double)(u_loc_single.map[k_first_scaling+i] * rightM0(type).operator()(i));
            }
        }
        else {
            long k_first_scaling= (k - numLeftBoundaryScalings - 1)*2 + inner_scaling_offset - 1;
            for (int i=innerM0.firstIndex(); i<=innerM0.lastIndex(); ++i) {
                u_loc_single_ptr=u_loc_single.map.find(k_first_scaling+i);
                if (u_loc_single_ptr!=u_loc_single_end) {
                    val += (*u_loc_single_ptr).second * innerM0(i);
                }
                //val += (long double)(u_loc_single.map[k_first_scaling+i] * innerM0(i));
            }
        }
    }
    return val;
}




template <typename PrimalBasis>
void
LocalRefinement<PrimalBasis>::test_reconstruct(int j, long k, XType xtype) const
{
    std::ofstream file("test.dat");
    if (xtype==XWavelet) {
        for (T x=0.; x<=1.; x+=0.001) {
            if (k<=numLeftBoundaryWavelets) {
                long k_first_scaling= basis.mra.rangeI(j+1).firstIndex() - 1;
                T val = 0.;
                for (int i=leftM1(k).firstIndex(); i<=leftM1(k).lastIndex(); ++i) {
                    val += basis.mra.phi(x,j+1,k_first_scaling + i,0) * leftM1(k).operator()(i);
                }
                file << x << " " << basis.psi(x,j,k,0) << " " << val << std::endl;
            }
            else if (k>basis.cardJ(j)-numRightBoundaryWavelets) {
                int type = basis.cardJ(j) - k + 1;
                long k_first_scaling= basis.mra.rangeI(j+1).lastIndex() - rightM1(type).length();
                T val = 0.;
                for (int i=rightM1(type).firstIndex(); i<=rightM1(type).lastIndex(); ++i) {
                    val += basis.mra.phi(x,j+1,k_first_scaling + i,0) * rightM1(type).operator()(i);
                }
                file << x << " " << basis.psi(x,j,k,0) << " " << val << std::endl;
            }
            else {
                long k_first_scaling= (k - numLeftBoundaryWavelets - 1)*2 + inner_wavelet_offset - 1;
                T factor = ((basis.d % 2 != 0) && (k>basis.cardJ(j)/2.)) ? -1. : 1.;
                if (x==0) {
                    std::cerr << "  -> " << k << " " << factor << std::endl;
                }
                T val = 0.;
                for (int i=innerM1.firstIndex(); i<=innerM1.lastIndex(); ++i) {
                    val += basis.mra.phi(x,j+1,k_first_scaling + i,0) * factor*innerM1(i);
                }
                file << x << " " << basis.psi(x,j,k,0) << " " <<  val << std::endl;
            }
        }
    }
    else {
        for (T x=0.; x<=1.; x+=0.001) {
            if (k<=basis.mra.rangeIL(j).lastIndex()) {
                int offset = basis.mra.rangeI(basis.j0).firstIndex()-1;
                long k_first_scaling= basis.mra.rangeI(j+1).firstIndex() - 1;
                T val = 0.;
                if (x==0) {
                    std::cerr << "Left boundary scaling." << std::endl;
                    std::cerr << k << " " << leftM0(k-offset) << std::endl;
                }
                for (int i=leftM0(k-offset).firstIndex(); i<=leftM0(k-offset).lastIndex(); ++i) {
                    val += basis.mra.phi(x,j+1,k_first_scaling + i,0) * leftM0(k-offset).operator()(i);
                }
                file << x << " " << basis.mra.phi(x,j,k,0) << " " << val << std::endl;
            }
            else if (k>=basis.mra.rangeIR(j).firstIndex()) {
                int type = basis.mra.rangeI(j).lastIndex() - k + 1;
                std::cerr << "    k = " << k << " " << type << std::endl;
                long k_first_scaling= basis.mra.rangeI(j+1).lastIndex() - rightM0(type).length();
                T val = 0.;
                if (x==0) {
                    std::cerr << "right boundary scaling." << std::endl;
                    std::cerr << k << " " << rightM0(type) << std::endl;
                }
                for (int i=rightM0(type).firstIndex(); i<=rightM0(type).lastIndex(); ++i) {
                    val += basis.mra.phi(x,j+1,k_first_scaling + i,0) * rightM0(type).operator()(i);
                }
                file << x << " " << basis.mra.phi(x,j,k,0) << " " << val << std::endl;
            }
            else {
                long k_first_scaling= (k - numLeftBoundaryScalings - 1)*2 + inner_scaling_offset - 1;
                T val = 0.;
                if (x==0) {
                    std::cerr << "Inner boundary scaling." << std::endl;
                    std::cerr << k << " " << innerM0 << std::endl;
                }
                for (int i=innerM0.firstIndex(); i<=innerM0.lastIndex(); ++i) {
                    val += basis.mra.phi(x,j+1,k_first_scaling + i,0) * innerM0(i);
                }
                file << x << " " << basis.mra.phi(x,j,k,0) << " " <<  val << std::endl;
            }
        }
    }
    file.close();
}

}   // namespace lawa
