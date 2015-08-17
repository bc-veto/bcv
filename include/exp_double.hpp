#ifndef OMEGA_EXP_DOUBLE_HPP
#define OMEGA_EXP_DOUBLE_HPP

#include <cmath>

#include "functors.hpp"

// a C99 function that exists in libm but is not declared in C89 <cmath>
double log1p(double);

namespace omega
{
    
    // we have to shadow std::log so we must provide a forwarding function
    
    const double log(const double& a)
    {
        assert(a >= 0);
        return std::log(a);
    }
    
    // a replacement for double, tailored for Bayesian marginalization
    // supports only positive numbers and stores their logarithm
    // achieves a vastly expanded range (e^(10^300))
    
    class exp_double
    {      
        
        double _; // internal storage
        
        // private constructor allowing friend functions to directly 
        // initialize the internal storage; a dummy argument is used to 
        // differentiate from conversion
        exp_double(const double& a, int) : _(a) {}    
        
    public:
        
        // default value is 0, internally represented as log(0) = -inf        
        exp_double() : _(std::log(0)) {}
        
        // default copy construction is correct
        // default destruction is correct
        // default assignment is correct
       
        // implicitly convert from double  
        exp_double(const double& a) : _(std::log(a)) 
        { 
            assert(a >= 0); 
        }
        
        // implicitly assign from double
        exp_double& operator=(const double& a)
        {
            _ = log(a);
            return *this;
        }
        
        // implicitly convert to double (DISABLED)
        //operator double() const
        //{
        //    assert(_ < 700);
        //    return std::exp(_);
        //}
        
        // allow omega:exp and omega::log access to the internals
        
        friend const exp_double exp(const double&);
        friend const double log(const exp_double&);
        
        // allow the computed assignments we can perform in-place to access
        // the internals
        
        friend exp_double& operator*=(exp_double&, const exp_double&);
        friend exp_double& operator*=(exp_double&, const double&);
        friend exp_double& operator/=(exp_double&, const exp_double&);
        friend exp_double& operator/=(exp_double&, const double&);
        
    }; // class exp_double
    
    // exp and log are closely integrated with the exp_double class
    
    // shadows std::exp(double)
    const exp_double exp(const double& a)
    {
        // uses friendship and dummy argument to access private constructor
        // that directly initializes the internal storage
        return exp_double(a, 0);
    }
    
    // overloads std::log
    const double log(const exp_double& a)
    {
        // uses friendship to access internal storage
        return a._;
    }
    
    // other functions use exp and log to interact with exp_double
    
    // output for debugging
    
    ostream& operator<<(ostream& a, const exp_double& b)
    {
        if ((-300 < log(b)) && (log(b) < 300))
            return a << std::exp(log(b));
        else
            return a << "exp(" << log(b) << ")";
    }
    
    // transcendental functions
    
    const exp_double sqrt(const exp_double& a)
    {
        return exp(log(a) / 2);
    }
        
    const exp_double pow(const exp_double& a, const double& b)
    {
        return exp(log(a) * b);
    }
    
    const exp_double pow(const exp_double& a, const int& b)
    {
        return exp(log(a) * b);
    }      
    
    // arithmetic operations
    
    // a * b = exp(log(a * b)) = exp(log(a) + log(b))
    
    const exp_double operator*(const exp_double& a, const exp_double& b)
    {
        return exp(log(a) + log(b));
    }
        
    const exp_double operator*(const exp_double& a, const double& b)
    {
        return exp(log(a) + log(b));
    }
    
    const exp_double operator*(const double& a, const exp_double& b)
    {
        return exp(log(a) + log(b));
    }
    
    // a / b = exp(log(a / b)) = exp(log(a) - log(b))
    
    const exp_double operator/(const exp_double& a, const exp_double& b)
    {
        return exp(log(a) - log(b));
    }
    
    const exp_double operator/(const exp_double& a, const double& b)
    {
        return exp(log(a) - log(b));
    }
        
    const exp_double operator/(const double& a, const exp_double& b)
    {
        return exp(log(a) - log(b));
    }   
    
    // addition is the most compex operation to implement;
    // to prevent overflow we compute via std::exp(-abs(a - b)), 
    // to improve accuracy we use log1p rather than log(1 + x)
    //
    // we use the identity
    // for 0 <= a <= b :
    // a + b = (a / b + 1) * b
    //       = exp(log((a / b + 1) * b))
    //       = exp(log(a / b + 1) + log(b))
    //       = exp(log(exp(log(a / b)) + 1) + log(b))
    //       = exp(log(exp(log(a) - log(b)) + 1) + log(b))
        
    const exp_double operator+(const exp_double& a, const exp_double& b)
    {
        if (log(a) < log(b))
        {
            return exp(log1p(std::exp(log(a) - log(b))) + log(b));
        }
        else if (log(b) < log(a))
        {
            return exp(log1p(std::exp(log(b) - log(a))) + log(a));
        }
        else
        {
            return exp(log(2) + log(a));
        }           
    }
    
    const exp_double operator+(const exp_double& a, const double& b)
    {
        return a + exp(log(b));
    }
    
    const exp_double operator+(const double& a, const exp_double& b)
    {
        return exp(log(a)) + b;
    }    
    
    // subtraction is simplified by requiring that a - b >= 0 ;
    // otherwise similar to addition
    
    const exp_double operator-(const exp_double& a, const exp_double& b)
    {
        assert(log(a) >= log(b));
        return exp(log1p(-std::exp(log(b) - log(a))) + log(a));
    }
    
    const exp_double operator-(const exp_double& a, const double& b)
    {
        return a - exp(log(b));
    }
    
    const exp_double operator-(const double& a, const exp_double& b)
    {
        return exp(log(a)) - b;
    }
    
    // computed assignment
    
    // computed assigned addition and subtraction cannot be implemented 
    // in-place, we forward to their respective binary operators
    
    exp_double& operator+=(exp_double& a, const exp_double& b)
    {
        return a = a + b;
    }

    exp_double& operator+=(exp_double& a, const double& b)
    {
        return a = a + b;
    }   
    
    exp_double& operator-=(exp_double& a, const exp_double& b)
    {
        return a = a - b;
    }

    exp_double& operator-=(exp_double& a, const double& b)
    {
        return a = a - b;
    }    
    
    // computed assigned multiplication and division can be implemented
    // in-place and are given access to the class internals
    
    exp_double& operator*=(exp_double& a, const exp_double& b)
    {        
        a._ += b._;
        return a;
    }
    
    exp_double& operator*=(exp_double& a, const double& b)
    {        
        assert(b >= 0);
        a._ += std::log(b);
        return a;
    }
    
    exp_double& operator/=(exp_double& a, const exp_double& b)
    {        
        a._ -= b._;
        return a;
    }
    
    exp_double& operator/=(exp_double& a, const double& b)
    {        
        assert(b > 0);
        a._ -= std::log(b);
        return a;
    }    
    
    // comparison operators
    //
    // when mixing with other types, we take into account the 
    // non-negativity of exp_doubles
    
    bool operator==(const exp_double& a, const exp_double& b)
    {
        return log(a) == log(b);
    }
    
    bool operator==(const double& a, const exp_double& b)
    {
        return (a >= 0) && (log(a) == log(b));
    }
   
    bool operator==(const exp_double& a, const double& b)
    {
        return (b >= 0) && (log(a) == log(b));
    }

    
    bool operator!=(const exp_double& a, const exp_double& b)
    {
        return log(a) != log(b);
    }
    
    bool operator!=(const double& a, const exp_double& b)
    {
        return (a < 0) || (log(a) != log(b));
    }
   
    bool operator!=(const exp_double& a, const double& b)
    {
        return (b < 0) || (log(a) != log(b));
    }
    
    
    bool operator<(const exp_double& a, const exp_double& b)
    {
        return log(a) < log(b);
    }
    
    bool operator<(const double& a, const exp_double& b)
    {
        return (a < 0) || (log(a) < log(b));
    }     
    
    bool operator<(const exp_double& a, const double& b)
    {
        return (b >= 0) && (log(a) < log(b));
    }
    
    
    bool operator<=(const exp_double& a, const exp_double& b)
    {
        return log(a) <= log(b);
    }
    
    bool operator<=(const double& a, const exp_double& b)
    {
        return (a <= 0) || (log(a) <= log(b));
    }     
    
    bool operator<=(const exp_double& a, const double& b)
    {
        return (b >= 0) && (log(a) <= log(b));
    }
    
    
    bool operator>(const exp_double& a, const exp_double& b)
    {
        return log(a) > log(b);
    }
    
    bool operator>(const double& a, const exp_double& b)
    {
        return (a > 0) && (log(a) > log(b));
    }     
    
    bool operator>(const exp_double& a, const double& b)
    {
        return (b < 0) || (log(a) > log(b));
    }   
    
    
    bool operator>=(const exp_double& a, const exp_double& b)
    {
        return log(a) >= log(b);
    }
    
    bool operator>=(const double& a, const exp_double& b)
    {
        return (a >= 0) && (log(a) >= log(b));
    }     
    
    bool operator>=(const exp_double& a, const double& b)
    {
        return (b <= 0) || (log(a) >= log(b));
    }    
    
    
    // tell the type inference system that doubles are promoted to 
    // exp_doubles in mixed-type arithmetic
    
    template<> class promoting_binary_functor<double, exp_double>
    {
    public:
        typedef exp_double type;
    };
    
    // tell the type inference system that exp(double) -> exp_double
        
    template<> class exp_functor<double>
    { 
    public: 
        typedef exp_double type; 
        const type operator()(const double& a) const { return exp(a); } 
    };
    
    // tell the type inference system that log(exp_double) -> double
    
    template<> class log_functor<exp_double>
    {
    public:
        typedef double type;
        const type operator()(const exp_double& a) const { return log(a); }
    };
    
} // namespace omega

#endif // OMEGA_EXP_DOUBLE_HPP
