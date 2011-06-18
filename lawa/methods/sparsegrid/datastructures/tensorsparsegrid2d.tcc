namespace lawa {

template <typename T, typename Basis2D, typename S1_x, typename S1_y, typename S2_x, typename S2_y>
TensorSparseGrid2D<T, Basis2D, S1_x, S1_y, S2_x, S2_y>::TensorSparseGrid2D
                                                        (const Basis2D &_basis,
                                                         const S1_x &_s1_x, const S1_y &_s1_y,
                                                         const S2_x &_s2_x, const S2_y &_s2_y,
                                                         int _J, T _gamma)
: basis(_basis), s1_x(_s1_x), s1_y(_s1_y), s2_x(_s2_x), s2_y(_s2_y),
  j0_x(basis.first.j0), j0_y(basis.second.j0), J(_J), gamma(_gamma),
  blockassembler1d(basis.first), dim(0)
{
    int fi=1, li=1;

    for (int ix=0; ix<=J; ++ix) {
        for (int iy=0; iy<=J; ++iy) {
            if (ix+iy>J) break;

            int *pair = new int[6];
            pair[0] = ix; pair[1] = iy;

            if (ix==0) {
                pair[4] = basis.first.mra.cardI(j0_x);
                if (iy==0) {
                    li = fi+basis.first.mra.cardI(j0_x)*basis.second.mra.cardI(j0_y)-1;
                    pair[5] = basis.second.mra.cardI(j0_y);
                }
                else {
                    li = fi+basis.first.mra.cardI(j0_x)*basis.second.cardJ(j0_y+iy-1)-1;
                    pair[5] = basis.second.cardJ(j0_y+iy-1);
                }
            }
            else {
                pair[4] = basis.first.cardJ(j0_x+ix-1);
                if (iy==0) {
                    li = fi+basis.first.cardJ(j0_x+ix-1)*basis.second.mra.cardI(j0_y)-1;
                    pair[5] = basis.second.mra.cardI(j0_y);
                }
                else {
                    li = fi+basis.first.cardJ(j0_x+ix-1)*basis.second.cardJ(j0_y+iy-1)-1;
                    pair[5] = basis.second.cardJ(j0_y+iy-1);
                }
            }
            pair[2] = fi; pair[3] = li;
            fi = li+1;
            dim += pair[4] * pair[5];
            sg_blocks.push_back(pair);
        }
    }

    int levelpair_pos = 0;
    for (int level1=0; level1<=J; ++level1) {
        for (int level2=0; level2<=J; ++level2) {
            levelpair_map[std::pair<int,int>(level1,level2)] = levelpair_pos;
            ++levelpair_pos;
        }
    }
    this->assembleMatrices();
}

template<typename T, typename Basis2D, typename S1_x, typename S1_y, typename S2_x, typename S2_y>
int
TensorSparseGrid2D<T,Basis2D,S1_x,S1_y,S2_x,S2_y>::getDimension() const
{
    return dim;
}

template<typename T, typename Basis2D, typename S1_x, typename S1_y, typename S2_x, typename S2_y>
IndexSet<Index2D>
TensorSparseGrid2D<T,Basis2D,S1_x,S1_y,S2_x,S2_y>::getIndexSet() const
{
    IndexSet<Index2D> sparsegridindexset;
    for (int block=0; block<(int)sg_blocks.size(); ++block) {
        int i1 = sg_blocks[block][0]-1;
        int i2 = sg_blocks[block][1]-1;

        if (i1<0) {
            for (int k1=basis.first.mra.rangeI(j0_x).firstIndex(); k1<=basis.first.mra.rangeI(j0_x).lastIndex(); ++k1) {
                if (i2<0) {
                    for (int k2=basis.second.mra.rangeI(j0_y).firstIndex(); k2<=basis.second.mra.rangeI(j0_y).lastIndex(); ++k2) {
                        Index2D index2d(Index1D(j0_x,k1,XBSpline),Index1D(j0_y,k2,XBSpline));
                        sparsegridindexset.insert(index2d);
                    }
                }
                else {
                    for (int k2=basis.second.rangeJ(j0_y+i2).firstIndex(); k2<=basis.second.rangeJ(j0_y+i2).lastIndex(); ++k2) {
                        Index2D index2d(Index1D(j0_x,k1,XBSpline),Index1D(j0_y+i2,k2,XWavelet));
                        sparsegridindexset.insert(index2d);
                    }
                }
            }
        }
        else {
            for (int k1=basis.first.rangeJ(j0_x+i1).firstIndex(); k1<=basis.first.rangeJ(j0_x+i1).lastIndex(); ++k1) {
                if (i2<0) {
                    for (int k2=basis.second.mra.rangeI(j0_y).firstIndex(); k2<=basis.second.mra.rangeI(j0_y).lastIndex(); ++k2) {
                        Index2D index2d(Index1D(j0_x+i1,k1,XWavelet),Index1D(j0_y,k2,XBSpline));
                        sparsegridindexset.insert(index2d);
                    }
                }
                else {
                    for (int k2=basis.second.rangeJ(j0_y+i2).firstIndex(); k2<=basis.second.rangeJ(j0_y+i2).lastIndex(); ++k2) {
                        Index2D index2d(Index1D(j0_x+i1,k1,XWavelet),Index1D(j0_y+i2,k2,XWavelet));
                        sparsegridindexset.insert(index2d);
                    }
                }
            }
        }
    }
    return sparsegridindexset;
}

template <typename T, typename Basis2D, typename S1_x, typename S1_y, typename S2_x, typename S2_y>
void
TensorSparseGrid2D<T, Basis2D, S1_x, S1_y, S2_x, S2_y>::toCoefficients
                                   (const DenseVectorT &vec,
                                    Coefficients<Lexicographical,T,Index2D> &sparsegridcoefficients)
{
    int pos=1;
    for (int block=0; block<(int)sg_blocks.size(); ++block) {
        int i1 = sg_blocks[block][0]-1;
        int i2 = sg_blocks[block][1]-1;

        //std::cout << "(" << sg_blocks[lp][0] << ", " << sg_blocks[lp][1] << ")" << std::endl;

        if (i1<0) {
            for (int k1=basis.first.mra.rangeI(j0_x).firstIndex(); k1<=basis.first.mra.rangeI(j0_x).lastIndex(); ++k1) {
                if (i2<0) {
                    for (int k2=basis.second.mra.rangeI(j0_y).firstIndex(); k2<=basis.second.mra.rangeI(j0_y).lastIndex(); ++k2) {
                        Index2D index2d(Index1D(j0_x,k1,XBSpline),Index1D(j0_y,k2,XBSpline));
                        sparsegridcoefficients[index2d] = vec(pos);
                        ++pos;
                    }
                }
                else {
                    for (int k2=basis.second.rangeJ(j0_y+i2).firstIndex(); k2<=basis.second.rangeJ(j0_y+i2).lastIndex(); ++k2) {
                        Index2D index2d(Index1D(j0_x,k1,XBSpline),Index1D(j0_y+i2,k2,XWavelet));
                        sparsegridcoefficients[index2d] = vec(pos);
                        ++pos;
                    }
                }
            }
        }
        else {
            for (int k1=basis.first.rangeJ(j0_x+i1).firstIndex(); k1<=basis.first.rangeJ(j0_x+i1).lastIndex(); ++k1) {
                if (i2<0) {
                    for (int k2=basis.second.mra.rangeI(j0_y).firstIndex(); k2<=basis.second.mra.rangeI(j0_y).lastIndex(); ++k2) {
                        Index2D index2d(Index1D(j0_x+i1,k1,XWavelet),Index1D(j0_y,k2,XBSpline));
                        sparsegridcoefficients[index2d] = vec(pos);
                        ++pos;
                    }
                }
                else {
                    for (int k2=basis.second.rangeJ(j0_y+i2).firstIndex(); k2<=basis.second.rangeJ(j0_y+i2).lastIndex(); ++k2) {
                        Index2D index2d(Index1D(j0_x+i1,k1,XWavelet),Index1D(j0_y+i2,k2,XWavelet));
                        sparsegridcoefficients[index2d] = vec(pos);
                        ++pos;
                    }
                }
            }
        }
    }
}

template<typename T, typename Basis2D, typename S1_x, typename S1_y, typename S2_x, typename S2_y>
flens::DenseVector<flens::Array<T> >
TensorSparseGrid2D<T,Basis2D,S1_x,S1_y,S2_x,S2_y>::operator*(const DenseVectorT &x) const
{
    std::vector<DenseMatrixT> y;
    for (int i=0; i<(int)sg_blocks.size(); ++i) {
        int dim_x = sg_blocks[i][4], dim_y = sg_blocks[i][5];
        DenseMatrixT Yi(dim_y,dim_x);
        y.push_back(Yi);
    }


    for (int col_block=0; col_block<(int)sg_blocks.size(); ++col_block) {
        int dim_x = sg_blocks[col_block][4], dim_y = sg_blocks[col_block][5];
        DenseMatrixT Xj(dim_y,dim_x);
        int fi = sg_blocks[col_block][2], li=fi+dim_y-1;
        for (int l=1; l<=dim_x; ++l) {
            Xj(_,l) = x(_(fi,li));
            fi = li+1;
            li = fi+dim_y-1;
        }

        for (int row_block=0; row_block<(int)sg_blocks.size(); ++row_block) {
            std::cout << "Multiplying (" << sg_blocks[row_block][0] << ", " << sg_blocks[row_block][1]
                      << ") * (" << sg_blocks[col_block][0] << ", " << sg_blocks[col_block][1] << ")" << std::endl;
            y[row_block] += this->block_multiplication(row_block,col_block,Xj);
        }
    }

    DenseVectorT ret(dim);
    int pos=1;
    for (int block=0; block<(int)sg_blocks.size(); ++block) {
        int dim_x = sg_blocks[block][4], dim_y = sg_blocks[block][5];
        for (int i=1; i<=dim_x; ++i) {
            for (int j=1; j<=dim_y; ++j) {
                ret(pos) = y[block](j,i);
                ++pos;
            }
        }
    }
    return ret;
}

template <typename T, typename Basis2D, typename S1_x, typename S1_y, typename S2_x, typename S2_y>
flens::GeMatrix<flens::FullStorage<T, cxxblas::ColMajor> >
TensorSparseGrid2D<T, Basis2D, S1_x, S1_y, S2_x, S2_y>::block_multiplication(int i, int j,
                                                                             const DenseMatrixT &Xj)
                                                                             const
{

    int dim_Lambda_x =       sg_blocks[i][4], dim_Lambda_y =       sg_blocks[i][5];
    int dim_Lambda_prime_x = sg_blocks[j][4], dim_Lambda_prime_y = sg_blocks[j][5];


    std::pair<int,int> lp_x(sg_blocks[i][0],sg_blocks[j][0]);
    std::pair<int,int> lp_y(sg_blocks[i][1],sg_blocks[j][1]);

    LevelPairMap::const_iterator it;
    it = levelpair_map.find(lp_x);
    int pos_x  = (*it).second;
    it = levelpair_map.find(lp_y);
    int pos_y  = (*it).second;

    int NumRows_x = matrixblocks_s1_x[pos_x].numRows();
    int NumCols_x = matrixblocks_s1_x[pos_x].numCols();
    int NumRows_y = matrixblocks_s1_y[pos_y].numRows();
    int NumCols_y = matrixblocks_s1_y[pos_y].numCols();

    assert(NumRows_x == dim_Lambda_x);
    assert(NumRows_y == dim_Lambda_y);
    assert(NumCols_x == dim_Lambda_prime_x);
    assert(NumCols_y == dim_Lambda_prime_y);

    assert(Xj.numCols() == dim_Lambda_prime_x);
    assert(Xj.numRows() == dim_Lambda_prime_y);


    DenseMatrixT Yi(dim_Lambda_y,dim_Lambda_x);

    if (dim_Lambda_y*dim_Lambda_prime_x <= dim_Lambda_prime_y*dim_Lambda_x) {
        std::cout << "->: " << dim_Lambda_y << " x " << dim_Lambda_prime_x << std::endl;
        DenseMatrixT S_x_otimes_I_y_X(dim_Lambda_y,dim_Lambda_prime_x);
        for (int j=1; j<=dim_Lambda_prime_x; ++j) {
            DenseVectorT tmp = matrixblocks_s1_y[pos_y]*Xj(_,j);
            S_x_otimes_I_y_X(_,j) = tmp;
        }

        for (int i=1; i<=dim_Lambda_y; ++i) {
            DenseVectorT v=S_x_otimes_I_y_X(i,_);
            DenseVectorT tmp = matrixblocks_s1_x[pos_x]*v;
            Yi(i,_) = tmp;
        }

        for (int j=1; j<=dim_Lambda_prime_x; ++j) {
            DenseVectorT tmp = matrixblocks_s2_y[pos_y]*Xj(_,j);
            S_x_otimes_I_y_X(_,j) = tmp;
        }

        for (int i=1; i<=dim_Lambda_y; ++i) {
            DenseVectorT v=S_x_otimes_I_y_X(i,_);
            DenseVectorT tmp = matrixblocks_s2_x[pos_x]*v;
            Yi(i,_) += tmp;
        }
    }

    else {
        std::cout << "->: " << dim_Lambda_prime_y << " x " << dim_Lambda_x << std::endl;
        DenseMatrixT I_x_otimes_S_y_X(dim_Lambda_prime_y,dim_Lambda_x);
        for (int i=1; i<=dim_Lambda_prime_y; ++i) {
            DenseVectorT v=Xj(i,_);
            DenseVectorT tmp = matrixblocks_s1_x[pos_x]*v;
            I_x_otimes_S_y_X(i,_) = tmp;
        }
        for (int j=1; j<=dim_Lambda_x; ++j) {
            DenseVectorT tmp = matrixblocks_s1_y[pos_y]*I_x_otimes_S_y_X(_,j);
            Yi(_,j) = tmp;
        }

        for (int i=1; i<=dim_Lambda_prime_y; ++i) {
            DenseVectorT v=Xj(i,_);
            DenseVectorT tmp = matrixblocks_s2_x[pos_x]*v;
            I_x_otimes_S_y_X(i,_) = tmp;
        }
        for (int j=1; j<=dim_Lambda_x; ++j) {
            DenseVectorT tmp = matrixblocks_s2_y[pos_y]*I_x_otimes_S_y_X(_,j);
            Yi(_,j) += tmp;
        }
    }

    return Yi;
}

template <typename T, typename Basis2D, typename S1_x, typename S1_y, typename S2_x, typename S2_y>
void
TensorSparseGrid2D<T, Basis2D, S1_x, S1_y, S2_x, S2_y>::assembleMatrices()
{
    for (LevelPairMap::const_iterator it=levelpair_map.begin(); it!=levelpair_map.end(); ++it) {
        int i1 = (*it).first.first;
        int i2 = (*it).first.second;
        int pos = (*it).second;
        //std::cout << "(" << i1 << ", " << i2 << "): " << pos << std::endl;
        SparseMatrixT A_s1_x = blockassembler1d.assembleBlock(s1_x,i1-1,i2-1,0.);
        SparseMatrixT A_s1_y = blockassembler1d.assembleBlock(s1_y,i1-1,i2-1,0.);
        SparseMatrixT A_s2_x = blockassembler1d.assembleBlock(s2_x,i1-1,i2-1,0.);
        SparseMatrixT A_s2_y = blockassembler1d.assembleBlock(s2_y,i1-1,i2-1,0.);
        matrixblocks_s1_x.push_back(A_s1_x);
        matrixblocks_s1_y.push_back(A_s1_y);
        matrixblocks_s2_x.push_back(A_s2_x);
        matrixblocks_s2_y.push_back(A_s2_y);
    }
}




}   // namespace lawa
