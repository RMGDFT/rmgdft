
#include "transition.h"

// Auxiliary function to get all prime factors of an integer
// less than a given value stop.
void GetPrimeFactors(std::vector<int>& factors, int val, int stop) {
   for(int i = 2;i <= stop;i++) {
       while( !(val % i) ) {
           factors.push_back(i);
           val = val / i;
       }
   }
}
