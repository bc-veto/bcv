#ifndef OMEGA_INTERVAL_HPP
#define OMEGA_INTERVAL_HPP

#include <cmath>
#include <limits>
#include <algorithm>

namespace omega
{
    
    const double pi = 3.14159265358979323846;
    
    using namespace std;
    
    using std::cos;
    using std::sin;
    using std::log;
    using std::exp;
    using ::log1p;
    using std::max;
    using std::min;
    
    class interval
    {
    public:
        double a, b;
        interval() : a(0), b(0) {};
        explicit interval(double a_) : a(a_), b(a_) {};
        interval(double a_, double b_) : a(a_), b(b_) {};
        interval& operator=(const double& a_) { a = b = a_; return *this; }
    };
    
    inline interval operator-(const interval& x)
    {
        return interval(-x.b, -x.a);
    }
    
    inline interval operator+(const interval& x, const interval& y)
    {
        return interval(x.a + y.a, x.b + y.b);
    }
    
    inline interval operator-(const interval& x, const interval& y)
    {
        return interval(x.a - y.b, x.b - y.a);
    }
    
    inline interval operator*(const interval& x, const interval& y)
    {
        double aa = x.a * y.a;
        double ab = x.a * y.b;
        double ba = x.b * y.a;
        double bb = x.b * y.b;        
        return interval(min(min(aa, ab), min(ba, bb)), max(max(aa, ab), max(ba, bb)));
    }
    
    inline interval inverse(const interval& x)
    {
        if ((x.a < 0) || (x.b > 0))
        {
            return interval(1./x.b, 1./x.a);
        }
        else
        {
            return interval(-std::numeric_limits<double>::infinity(), std::numeric_limits<double>::infinity());
        }
    }
    
    inline interval operator/(const interval& x, const interval& y)
    {
        return x * inverse(y);
    }
    
    inline interval cos(const interval& x)
    {
        
        if (x.b < 0)
            return cos(-x);
        
        if (x.a < 0)
            return cos(interval(0, max(-x.a, x.b)));
        
        if (x.a >= 2 * pi)
        {
            double a = fmod(x.a, 2 * pi);
            return cos(interval(a, a + (x.b - x.a)));
        }
                        
        if (x.a >= pi)
            return -cos(x - interval(pi));
        
        if (x.a >= pi/2)
            return -cos(interval(pi) - x);
        
        if (x.b >= 2 * pi)
            return interval(-1, 1);
        
        if (x.b >= pi)
            return interval(-1, std::max(std::cos(x.a), std::cos(x.b)));
        
        return interval(std::cos(x.b), std::cos(x.a));
        
    }
    
    inline interval sin(const interval& x)
    {
        return cos(x - interval(pi/2));
    }
    
    inline interval log(const interval& x)
    {
        return interval(log(x.a), log(x.b));
    }
    
    inline interval exp(const interval& x)
    {
        return interval(exp(x.a), exp(x.b));
    }
    
    inline interval log1p(const interval& x)
    {
        return interval(log1p(x.a), log1p(x.b));
    }
    
    inline interval max(const interval& x, const interval& y)
    {
        return interval(max(x.a, y.a), max(x.b, y.b));
    }
    
    inline interval min(const interval& x, const interval& y)
    {
        return interval(min(x.a, y.a), min(x.b, y.b));
    }
    
    inline interval acos(const interval& x)
    {
        return interval(x.b, x.a);
    }
    
    inline double midpoint(const interval& x)
    {
        return (x.a + x.b) / 2;
    }
    
    inline double ith(const interval& x, int i, int n)
    {
        return (x.a * (i - n) + x.b * i) / n;
    }
    
    inline interval partition(const interval& x, int i, int n)
    {
        return interval((x.a * (n - i) + x.b * i) / n, (x.a * (n - i - 1) + x.b * (i + 1)) / n);
    }
    
    inline double length(const interval& x)
    {
        return x.b - x.a;
    }
    
    inline interval& operator+=(interval& x, const interval& y)
    {
        x.a += y.a;
        x.b += y.b;
        return x;
    }
    
    inline interval hull(const interval& x, const interval& y)
    {
        return interval(min(x.a, y.a), max(x.b, y.b));
    }
    
    inline interval hull(const interval& x, const double& y)
    {
        return hull(x, interval(y));
    }
    
    inline interval& operator*=(interval& x, const double& y)
    {
        if (y >= 0)
        {
            x.a *= y;
            x.b *= y;
        }
        else
        {
            x = interval(x.b * y, x.a * y);
        }
        return x;
    }
    
};

#endif // OMEGA_INTERVAL_HPP
