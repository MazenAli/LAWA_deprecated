#include <iostream>
#include <lawa/lawa.h>
#include <lawa/periodic/periodicsupport.h>

using namespace lawa;
using namespace std;

typedef double T;

int 
main(int arc, char* argv[]){
    
    PeriodicSupport<T> supp1(0.2, 0.8);
    PeriodicSupport<T> supp1b(0.25, 0.75);
    PeriodicSupport<T> supp2(0,1, 0.3, 0.7);
    PeriodicSupport<T> supp3(0,1, 0.2, 0.5);
    PeriodicSupport<T> supp4(0,1, 0.2, 0.9);
    PeriodicSupport<T> supp5(0,1, 0.4, 0.9);
    
    
    /* ======= OVERLAPS =================*/
    
    Support<T> common;
    
    cout << "Supports: " << supp1 << " and " << supp2 << endl;
    overlap(supp1, supp2, common);
    cout << "Overlap = " << overlap(supp1, supp2) << " : " << common << endl;
    cout << "Distance: " << lawa::distance(supp1, supp2) << endl;
    cout << "Minimal Overlap = " << minimal_overlap(supp1, supp2) << endl << endl;
    
    cout << "Supports: " << supp1 << " and " << supp3 << endl;
    overlap(supp1, supp3, common);
    cout << "Overlap = " << overlap(supp1, supp3) << " : " << common << endl;
    cout << "Distance: " << lawa::distance(supp1, supp3) << endl;
    cout << "Minimal Overlap = " << minimal_overlap(supp1, supp3) << endl << endl;
    
    cout << "Supports: " << supp1 << " and " << supp4 << endl;
    overlap(supp1, supp4, common);
    cout << "Overlap = " << overlap(supp1, supp4) << " : " << common << endl;
    cout << "Distance: " << lawa::distance(supp1, supp4) << endl;
    cout << "Minimal Overlap = " << minimal_overlap(supp1, supp4) << endl << endl;
    
    cout << "Supports: " << supp1 << " and " << supp5 << endl;
    overlap(supp1, supp5, common);
    cout << "Overlap = " << overlap(supp1, supp5) << " : " << common << endl;
    cout << "Distance: " << lawa::distance(supp1, supp5) << endl;
    cout << "Minimal Overlap = " << minimal_overlap(supp1, supp5) << endl << endl;
    
    cout << "Supports: " << supp2 << " and " << supp1b << endl;
    overlap(supp2, supp1b, common);
    cout << "Overlap = " << overlap(supp2, supp1b) << " : " << common << endl;
    cout << "Distance: " << lawa::distance(supp2, supp1b) << endl;
    cout << "Minimal Overlap = " << minimal_overlap(supp2, supp1b) << endl << endl;
    
    cout << "Supports: " << supp3 << " and " << supp1b << endl;
    overlap(supp3, supp1b, common);
    cout << "Overlap = " << overlap(supp3, supp1b) << " : " << common << endl;
    cout << "Distance: " << lawa::distance(supp3, supp1b) << endl;
    cout << "Minimal Overlap = " << minimal_overlap(supp3, supp1b) << endl << endl;
    
    cout << "Supports: " << supp4 << " and " << supp1b << endl;
    overlap(supp4, supp1b, common);
    cout << "Overlap = " << overlap(supp4, supp1b) << " : " << common << endl;
    cout << "Distance: " << lawa::distance(supp4, supp1b) << endl;
    cout << "Minimal Overlap = " << minimal_overlap(supp4, supp1b) << endl << endl;
    
    cout << "Supports: " << supp5 << " and " << supp1b << endl;
    overlap(supp5, supp1b, common);
    cout << "Overlap = " << overlap(supp5, supp1b) << " : " << common << endl;
    cout << "Distance: " << lawa::distance(supp5, supp1b) << endl;
    cout << "Minimal Overlap = " << minimal_overlap(supp5, supp1b) << endl << endl;
    
    cout << "Supports: " << supp2 << " and " << supp3 << endl;
    overlap(supp2, supp3, common);
    cout << "Overlap = " << overlap(supp2, supp3) << " : " << common << endl;
    cout << "Distance: " << lawa::distance(supp2, supp3) << endl;
    cout << "Minimal Overlap = " << minimal_overlap(supp2, supp3) << endl << endl;
    
    cout << "Supports: " << supp2 << " and " << supp4 << endl;
    overlap(supp2, supp4, common);
    cout << "Overlap = " << overlap(supp2, supp4) << " : " << common << endl;
    cout << "Distance: " << lawa::distance(supp2, supp4) << endl;
    cout << "Minimal Overlap = " << minimal_overlap(supp2, supp4) << endl << endl;
    
    cout << "Supports: " << supp2 << " and " << supp5 << endl;
    overlap(supp2, supp5, common);
    cout << "Overlap = " << overlap(supp2, supp5) << " : " << common << endl;
    cout << "Distance: " << lawa::distance(supp2, supp5) << endl;
    cout << "Minimal Overlap = " << minimal_overlap(supp2, supp5) << endl << endl;
    
    /* ======= SHIFTS & FACTORS =================*/
 
    T shift1 = 3;
    T shift2 = 0.75;
    T shift3 = 0.15;
        
    cout << "Support: " << supp1 << endl;
    cout << "   + " << shift1 << " : " << supp1 + shift1 << endl;
    cout << "   + " << shift2 << " : " << supp1 + shift2 << endl;
    cout << "   + " << shift3 << " : " << supp1 + shift3 << endl;
    
    cout << "Support: " << supp2 << endl;
    cout << "   + " << shift1 << " : " << supp2 + shift1 << endl;
    cout << "   + " << shift2 << " : " << supp2 + shift2 << endl;
    cout << "   + " << shift3 << " : " << supp2 + shift3 << endl;

    cout << "Support: " << supp3 << endl;
    cout << "   + " << shift1 << " : " << supp3 + shift1 << endl;
    cout << "   + " << shift2 << " : " << supp3 + shift2 << endl;
    cout << "   + " << shift3 << " : " << supp3 + shift3 << endl;
  
    cout << "Support: " << supp4 << endl;
    cout << "   + " << shift1 << " : " << supp4 + shift1 << endl;
    cout << "   + " << shift2 << " : " << supp4 + shift2 << endl;
    cout << "   + " << shift3 << " : " << supp4 + shift3 << endl;
    
    cout << "Support: " << supp5 << endl;
    cout << "   + " << shift1 << " : " << supp5 + shift1 << endl;
    cout << "   + " << shift2 << " : " << supp5 + shift2 << endl;
    cout << "   + " << shift3 << " : " << supp5 + shift3 << endl;
    
    
    PeriodicSupport<T> supp6(0.2, 0.6, 0.1, 0.5);
    
    return 0;
}