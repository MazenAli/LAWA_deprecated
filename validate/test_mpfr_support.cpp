#include <boost/math/bindings/mpfr.hpp>
#include <boost/math/special_functions/gamma.hpp>
#include <boost/math/special_functions/binomial.hpp>
#include <boost/math/special_functions/factorials.hpp>

using boost::math::policies::policy;
using boost::math::policies::digits10;

int main()
{
    typedef boost::math::policies::policy<boost::math::policies::digits2<64> > my_prec_policy;

   mpfr_class::set_dprec(500); // 500 bit precision
   //
   // Note that the argument to tgamma is an expression template,
   // that's just fine here:
   //
   /*
   mpfr_class v1 = boost::math::tgamma(sqrt(mpfr_class(2)));
   std::cout << std::setprecision(50) << v1 << std::endl;

   std::cout << std::endl;
   */

   std::cout << std::endl;

   //z.set_prec(256), x.set_prec(256), v2.set_prec(256);

   mpfr_class test = sqrt(mpfr_class(2));

   double test_d = mpfr_get_d(test.get_mpfr_t(), GMP_RNDN);

   std::cout << "1.4142135623730950488016887242096980785696718753769" << std::endl;
   std::cout << std::setprecision(50) << sqrt(mpfr_class(2)) << std::endl;
   std::cout << std::setprecision(50) << test_d << std::endl;
   std::cout << std::setprecision(50) << std::sqrt(2.L) << std::endl << std::endl << std::endl;

   std::cout << "0.88658142871925912508091761239199431116828551763805" << std::endl;
   std::cout << std::setprecision(50) << boost::math::tgamma(sqrt(mpfr_class(2)),policy<digits10<20> >()) << std::endl;
   std::cout << std::setprecision(50) << boost::math::tgamma((long double)sqrt(2.L)) << std::endl << std::endl;

   mpfr_class z,x;
   z = 0.25;
   x = 0.125;
   //std::cout << std::setprecision(50) << x << std::endl;
   //std::cout << std::setprecision(50) << z << std::endl;
   std::cout << "1.3046495990274244829148742432632630302142654055017" << std::endl;
   std::cout << boost::math::tgamma(z,x) << std::endl;
   std::cout << boost::math::tgamma((long double)0.25,(long double)0.125) << std::endl;

}
