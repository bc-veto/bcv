#include <algorithm>
#include <complex>
#include <iostream>
#include <limits>
#include <list>

// Omega

#include "exp_double.hpp"

namespace omega
{
    /*
    bool approximates(double a, double b)
    {
        return abs(a - b) < (std::max(a, b) * 1e-7);  
    }
     */
    
    template<typename A, typename B> bool approximates(A a,  B b)
    {
        return abs(log(a) - log(b)) < 1e-7;
    }    
    
    int failures(0);
    
    void testf(bool a, const char* b)
    {
        if (a)
        {
            cout << "pass: " << b << endl;
        }
        else
        {
            cout << "fail: " << b << endl;
            ++failures;
        }
    }
    
    #define TEST(A) cout << ((A) ? "pass: " : "fail: ") << #A << endl;
    
    #define SHOW(A) cout << "      " << #A << endl; A;
    
    void test()
    {

        TEST(approximates(log(2), std::log(2)))
        TEST(approximates(exp_double(1), 1))
        TEST(approximates(exp(0), 1))
        TEST(approximates(log(exp(1)), 1))
        TEST(approximates(sqrt(exp_double(4)), 2))
        TEST(approximates(pow(exp_double(2), 2), 4))
        TEST(approximates(pow(exp_double(2), 2.0), 4))
        
        TEST(approximates(exp_double(1) + 1, 2))
        TEST(approximates(1 + exp_double(1), 2))
        TEST(approximates(exp_double(1) + exp_double(1), 2))
        
        TEST(approximates(exp_double(2) - 1, 1))
        TEST(approximates(2 - exp_double(1), 1))
        TEST(approximates(exp_double(2) - exp_double(1), 1))

        TEST(approximates(exp_double(2) * 3, 6))
        TEST(approximates(2 * exp_double(3), 6))
        TEST(approximates(exp_double(2) * exp_double(3), 6))

        TEST(approximates(exp_double(6) / 3, 2))
        TEST(approximates(6 / exp_double(3), 2))
        TEST(approximates(exp_double(6) / exp_double(3), 2))
       
        SHOW(exp_double a(1))

        TEST(approximates(a += 1, 2))
        TEST(approximates(a -= 1, 1))
        TEST(approximates(a += exp_double(1), 2))
        TEST(approximates(a -= exp_double(1), 1))
        TEST(approximates(a *= 2, 2))
        TEST(approximates(a /= 2, 1))
        TEST(approximates(a *= exp_double(2), 2))
        TEST(approximates(a /= exp_double(2), 1))
        
        TEST(  exp_double(1) == exp_double(1))
        TEST(!(exp_double(1) == exp_double(2)))
        TEST(  (1) == exp_double(1))
        TEST(!((1) == exp_double(2)))
        TEST(  exp_double(1) == (1))
        TEST(!(exp_double(1) == (2)))
        
        TEST(  exp_double(1) != exp_double(2))
        TEST(!(exp_double(1) != exp_double(1)))
        TEST(  (1) != exp_double(2))
        TEST(!((1) != exp_double(1)))
        TEST(  exp_double(1) != (2))
        TEST(!(exp_double(1) != (1)))
        
        TEST(  exp_double(1) < exp_double(2))
        TEST(!(exp_double(2) < exp_double(1)))
        TEST(!(exp_double(2) < exp_double(2)))
        TEST(  (1) < exp_double(2))
        TEST(!((2) < exp_double(1)))
        TEST(!((2) < exp_double(2)))
        TEST(  exp_double(1) < (2))
        TEST(!(exp_double(2) < (1)))
        TEST(!(exp_double(2) < (2)))
       
        TEST(  exp_double(1) <= exp_double(2))
        TEST(!(exp_double(2) <= exp_double(1)))
        TEST(  exp_double(2) <= exp_double(2))    
        TEST(  (1) <= exp_double(2))
        TEST(!((2) <= exp_double(1)))
        TEST(  (2) <= exp_double(2))    
        TEST(  exp_double(1) <= (2))
        TEST(!(exp_double(2) <= (1)))
        TEST(  exp_double(2) <= (2))  
        
        TEST(  exp_double(2) > exp_double(1))
        TEST(!(exp_double(1) > exp_double(2)))
        TEST(!(exp_double(2) > exp_double(2)))
        TEST(  (2) > exp_double(1))
        TEST(!((1) > exp_double(2)))
        TEST(!((2) > exp_double(2)))
        TEST(  exp_double(2) > (1))
        TEST(!(exp_double(1) > (2)))
        TEST(!(exp_double(2) > (2)))
        
        TEST(  exp_double(2) >= exp_double(1))
        TEST(!(exp_double(1) >= exp_double(2)))
        TEST(  exp_double(2) >= exp_double(2))    
        TEST(  (2) >= exp_double(1))
        TEST(!((1) >= exp_double(2)))
        TEST(  (2) >= exp_double(2))    
        TEST(  exp_double(2) >= (1))
        TEST(!(exp_double(1) >= (2)))
        TEST(  exp_double(2) >= (2))   
        
        cout << endl << failures << " failures" << endl;
    }
}

using namespace omega;
using namespace std;

int main(int argc, char** argv)
{
   omega::test();
}
