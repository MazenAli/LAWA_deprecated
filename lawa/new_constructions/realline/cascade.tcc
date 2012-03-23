#ifndef LAWA_CONSTRUCTIONS_REALLINE_CASCADE_TCC
#define LAWA_CONSTRUCTIONS_REALLINE_CASCADE_TCC 1

//TODO: complete review of implementation style.
#include <cstdlib>
#include <cstdio>
#include <cmath>

#include <lawa/math/math.h>

namespace lawa {

template <typename T>
struct cascade
{
  T* vector;
  int level;
  int start;
  int end;
};

template <typename T>
T
cascade_ref(cascade<T>* c, int n)
{
  int h, i;

  h = pow2i<T>(c->level);
  if(n<h*c->start)
    return 0;
  if(n>h*c->end)
    return 0;
  i=n-h*c->start;
  return (c->vector)[i];
}

template <typename T>
void
cascade_output(cascade<T>* c)
{
  int i;

  for(i=c->start*pow2i<T>(c->level); i<=c->end*pow2i<T>(c->level); i++)
    {
      printf("%f ", cascade_ref(c,i));
    }
  printf("%i %i %i\n\n",c->level,c->start,c->end );
}

template <typename T>
cascade<T>*
eval(const T* mask, int /*length*/, int start, int end, int level)
{
  cascade<T>* old, *newv;
  int i,k,m,n;
  T sum,h;
  T* v;

  old=(cascade<T>*)malloc(sizeof(cascade<T>));
  newv=(cascade<T>*)malloc(sizeof(cascade<T>));
  old->level=0;
  old->start=0;
  old->end=0;
  old->vector=(T*)malloc(1*sizeof(T));
  std::uninitialized_fill_n(old->vector,1,T(0));  
  (old->vector)[0]=1;

  for(m=1; m<=level; m++)
    {
      newv->level=m;
      newv->start=start;
      newv->end=end;
      newv->vector=(T*)malloc((1+(end-start)*pow2i<T>(m))*sizeof(T));
      std::uninitialized_fill_n(newv->vector,1+(end-start)*pow2i<T>(m),T(0));

      for(n=start*pow2i<T>(m); n<=end*pow2i<T>(m); n++)
    {
      sum=0;
      for(k=start; k<=end; k++)
        {
          h=mask[k-start]*
        cascade_ref(old,n-k*pow2i<T>(m-1));
        

          sum+=h;
        }
      (newv->vector)[n-pow2i<T>(m)*start]=sum;
    }
      if (m>=2) {
          for(i=0; i<1+(end-start)*pow2i<T>(m-1);i++) {
              (old->vector)[i].~T();
          }
          free(old->vector);
      }
      old=(cascade<T>*)malloc(sizeof(cascade<T>));
      old->level=m;
      old->start=start;
      old->end=end;
      v=(T*)malloc((1+(end-start)*pow2i<T>(m))*sizeof(T));
      std::uninitialized_fill_n(v,1+(end-start)*pow2i<T>(m),T(0));
      for(i=0; i<1+(end-start)*pow2i<T>(m);i++) {
          v[i]=(newv->vector)[i];
          (newv->vector)[i].~T();
      }
      free(newv->vector);
      old->vector=v;
    }
    return old;
}

template <typename X, typename Y>
void
evalAtDyadicGrid_Cascade(const DenseVector<X> &sf, int J,
                         DenseVector<Y> &scaling)
{
    typedef typename X::ElementType T;
    
    cascade<T> *c = eval(sf.engine().data(), 0, sf.firstIndex(), sf.lastIndex(), J);

    scaling.resize((pow2i<T>(J))*(sf.lastIndex()-sf.firstIndex())+1, 0, T(0));
    for (int i=scaling.firstIndex(); i<=scaling.lastIndex(); ++i) {
        scaling(i) = c->vector[i];
    }
    scaling.engine().changeIndexBase(pow2i<T>(J)*sf.firstIndex());
}

} // namespace lawa

#endif // LAWA_CONSTRUCTIONS_REALLINE_CASCADE_TCC
