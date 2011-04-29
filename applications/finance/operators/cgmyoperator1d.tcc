namespace lawa {

template <typename T, typename Basis1D>
FinanceOperator1D<T, CGMY, Basis1D>::FinanceOperator1D(const Basis1D& _basis,
                                                     const Parameters<T,CGMY> &_params,
                                                     T _convection, T _reaction, T _R1, T _R2)
    : basis(_basis), params(_params), convection(_convection), reaction(_reaction),
      kernel(params), R1(_R1), R2(_R2),
      OneDivSqrtR2pR1(1./std::sqrt(R2+R1)), OneDivR2pR1(1./(R2+R1)), R1DivR1pR2(R1/(R1+R2)),
      integral(basis,basis)
{

}

template <typename T, typename Basis1D>
T
FinanceOperator1D<T, CGMY, Basis1D>::operator()(XType xtype1, int j1, int k1,
                                              XType xtype2, int j2, int k2) const
{
   typedef typename std::map<T,T>::const_iterator const_it;

   T int_val = 0.;

   GeMatrix<FullStorage<T,ColMajor> > varphi_row_deltas, varphi_col_deltas;
   varphi_row_deltas = computeDeltas<T,Basis1D>(basis,j1,k1,xtype1);
   varphi_col_deltas = computeDeltas<T,Basis1D>(basis,j2,k2,xtype2);

   if (R1!=0 && R2!=1) {
       varphi_row_deltas(_,1)  *=(R1+R2);  varphi_row_deltas(_,1)-=R1;
       varphi_row_deltas(_,2)  *= OneDivR2pR1*OneDivSqrtR2pR1;
       varphi_col_deltas(_,1)  *=(R1+R2);  varphi_col_deltas(_,1)-=R1;
       varphi_col_deltas(_,2)  *= OneDivR2pR1*OneDivSqrtR2pR1;
   }

   std::cout << varphi_row_deltas << std::endl;
   std::cout << varphi_col_deltas << std::endl;

   T part1=0., part2=0.;
   if (basis.d==2) {
       for (int lambda=varphi_row_deltas.rows().firstIndex();
                lambda<=varphi_row_deltas.rows().lastIndex(); ++lambda) {

           T x = varphi_row_deltas(lambda,1);
           if (fabs(varphi_row_deltas(lambda,2)) < 1e-14) continue;

           part1 += OneDivSqrtR2pR1*varphi_row_deltas(lambda,2)
                    *basis.generator(xtype2)((x+R1)* OneDivR2pR1,j2,k2,0)*kernel.c3;


           for (int mu=varphi_col_deltas.rows().firstIndex();
                    mu<=varphi_col_deltas.rows().lastIndex(); ++mu) {
               if (fabs(varphi_col_deltas(mu,2)) < 1e-14) continue;

               T y = varphi_col_deltas(mu,1);
               T c = varphi_col_deltas(mu,2)*varphi_row_deltas(lambda,2);


               if (fabs(x-y)>1e-10)  {
                   T value_tailintegral;
                   const_it it   = values_tailintegral.find(y-x);
                   const_it last = values_tailintegral.end();
                   if (it != last) value_tailintegral = (*it).second;
                   else {
                       value_tailintegral=kernel.ForthTailIntegral(y-x);
                       values_tailintegral[y-x] = value_tailintegral;
                   }
                   if (y-x>0)  part2 += c * (value_tailintegral - kernel.constants[2]);
                   else        part2 += c * (value_tailintegral - kernel.constants[3]);
               }
           }
       }
       int_val = part1 + part2;
   }

   else if (basis.d==3) {

       bool is_realline_basis=flens::IsSame<Basis1D, Basis<T,Primal,R,CDF> >::value;
       assert(is_realline_basis);

       int_val -= kernel.c3*integral(j1,k1,xtype1,1,j2,k2,xtype2,1);

       for (int mu=varphi_col_deltas.rows().firstIndex();
                mu<=varphi_col_deltas.rows().lastIndex(); ++mu) {

           T y = varphi_col_deltas(mu,1);
           if (fabs(varphi_col_deltas(mu,2)) < 1e-14) continue;

           part1 += varphi_col_deltas(mu,2)*basis.generator(xtype1)(y,j1,k1,0)*kernel.c4;
           part1 -= varphi_col_deltas(mu,2)*basis.generator(xtype1)(y,j1,k1,0)*kernel.c5;

           for (int lambda=varphi_row_deltas.rows().firstIndex();
                    lambda<=varphi_row_deltas.rows().lastIndex(); ++lambda) {

               if (fabs(varphi_row_deltas(lambda,2)) < 1e-14) continue;
               T x = varphi_row_deltas(lambda,1);
               T c = varphi_col_deltas(mu,2)*varphi_row_deltas(lambda,2);
               if (x!=y)  {
                   if (y-x>0)  {
                       part2 -= c * (kernel.SixthTailIntegral(y-x) - kernel.constants[6]);
                   }
                   else    {
                       part2 -= c * (kernel.SixthTailIntegral(y-x) - kernel.constants[7]);
                   }
               }
           }
       }
       int_val += part1 + part2;
   }
   else {
       assert(0);
   }
   return int_val;
}

}    //namespace lawa
