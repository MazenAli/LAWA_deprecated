#include <boost/math/bindings/mpfr.hpp>
#include <boost/math/special_functions/gamma.hpp>

int main()
{
   mpfr_class::set_dprec(500); // 500 bit precision
   //
   // Note that the argument to tgamma is an expression template,
   // that's just fine here:
   //
   mpfr_class v = boost::math::tgamma(sqrt(mpfr_class(2)));
   std::cout << std::setprecision(50) << v << std::endl;
   std::cout << std::setprecision(50) << boost::math::tgamma(sqrt(2.)) << std::endl;
   std::cout << std::setprecision(50) << boost::math::tgamma((long double)sqrt(2.L)) << std::endl;
}
