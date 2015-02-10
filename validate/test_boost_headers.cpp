#include <iostream>

#include <boost/math/special_functions/gamma.hpp>

#ifdef TRONE
    #include <tr1/unordered_map>
#elif BOOST
    #include <boost/unordered_set.hpp>
#else
    #include <ext/hash_set>
#endif
//#include <boost/unordered_map.hpp>


//using namespace boost::math;
using namespace std;
//using namespace lawa;

typedef double T;

int main (int argc, char *argv[]) {

    std::cerr << "Only testing for include errors due to boost." << std::endl;
    return 0;
}
