#ifndef OMEGA_FUNCTORS_HPP
#define OMEGA_FUNCTORS_HPP

#include <cassert>
#include <cmath>

#include <complex>

namespace omega
{
    
    using namespace std;
    
    using std::abs;
    
    template<class T> class abs_functor 
    { 
    public: 
        typedef T type; 
        const T operator()(const T& a) const { return abs(a); } 
    };
    
    template<class T> class abs_functor<complex<T> >
    { 
    public: 
        typedef T type; 
        const T operator()(const complex<T>& a) const { return abs(a); } 
    };    
    
    template<class T> class exp_functor 
    { 
    public: 
        typedef T type; 
        const T operator()(const T& a) const { return exp(a); } 
    };
    
    template<class T> class log_functor 
    { 
    public: 
        typedef T type; 
        const T operator()(const T& a) const { return log(a); } 
    };
    
    #define UNARY_OPERATION(NAME) \
    using std::NAME; \
    template<class T> class NAME##_functor \
    { \
    public: \
        typedef T type; \
        const T operator()(const T& a) const { return NAME(a); } \
    }; \
    
    UNARY_OPERATION(acos)
    UNARY_OPERATION(asin)
    UNARY_OPERATION(atan)
    UNARY_OPERATION(cos)
    UNARY_OPERATION(conj)
    UNARY_OPERATION(cosh)
    UNARY_OPERATION(log10)
    UNARY_OPERATION(sin)
    UNARY_OPERATION(sinh)
    UNARY_OPERATION(sqrt)
    UNARY_OPERATION(tan)
    UNARY_OPERATION(tanh) 
    
    #undef UNARY_OPERATION   

    #define COMPLEX_UNARY_OPERATION(NAME) \
    using std::NAME; \
    template<class T> class NAME##_functor; \
    template<class T> class NAME##_functor<complex<T> > \
    { \
    public: \
        typedef T type; \
        const T operator()(const complex<T>& a) const { return NAME(a); } \
    };
    
    COMPLEX_UNARY_OPERATION(real)
    COMPLEX_UNARY_OPERATION(imag)
    COMPLEX_UNARY_OPERATION(arg)
    COMPLEX_UNARY_OPERATION(norm)   
    
    #undef COMPLEX_UNARY_OPERATION    
    
    // promoting_binary_functor<T, U>::type is the return type of A (binary operator) B
    template<class T, class U> class promoting_binary_functor
    { 
        public:
            typedef T type; 
    };
            
    // if one of the types is real and the other complex, the result is complex
    template<class T, class U> class promoting_binary_functor<complex<T>, U>
    {
    public:
        typedef complex<typename promoting_binary_functor<T, U>::type> type;
    };
     
    template<class T, class U> class promoting_binary_functor<T, complex<U> >
    {
    public:
        typedef complex<typename promoting_binary_functor<T, U>::type> type;
    };      
    
    #define BINARY_OPERATION(NAME, OPERATOR)\
    \
    template<class T, class U> class NAME##_functor : public promoting_binary_functor<T, U> \
    { \
    public:\
        const typename promoting_binary_functor<T, U>::type operator()(const T& a, const U& b) const { return a OPERATOR b; } \
    }; \
        
    BINARY_OPERATION(plus, +)
    BINARY_OPERATION(minus, -)
    BINARY_OPERATION(times, *)
    BINARY_OPERATION(divided_by, /)
    
    #undef BINARY_OPERATION    
    
    #define BINARY_FUNCTION(NAME)\
    \
    using std:: NAME ; \
    template<class T, class U> class NAME##_functor : public promoting_binary_functor<T, U> \
    { \
    public:\
        const typename promoting_binary_functor<T, U>::type operator()(const T& a, const U& b) const { return NAME(a, b); } \
    }; \
    
    BINARY_FUNCTION(pow)
    BINARY_FUNCTION(atan2)
    
    #undef BINARY_FUNCTION    
    
    template<class U> class pow_functor<int, U> : public promoting_binary_functor<double, U> 
    { 
    public:
        const typename promoting_binary_functor<double, U>::type operator()(const int& a, const U& b) const { return pow(a, b); } \
    };
    
    template<class F, class T, class U> class bind1st
    {
        const F f_;
        const T& a_;
    public:
        typedef typename F::type type;
        bind1st(const F& f, const T& a) : f_(f), a_(a) {}
        const type operator()(const U& b) const { return f_(a_, b); }
    };
    
    template<class F, class T, class U> class bind2nd
    {
        const F f_;
        const U& b_;
    public:
        typedef typename F::type type;
        bind2nd(const F& f, const U& b) : f_(f), b_(b) {}
        const type operator()(const T& a) const { return f_(a, b_); }        
    };
    
    // supplement the complex operators
   
    template<typename T> const complex<T> operator+(const complex<T>& a, const int b) { return a + T(b); } 
    template<typename T> const complex<T> operator+(const int a, const complex<T>& b) { return T(a) + b; }  
    template<typename T> const complex<T> operator-(const complex<T>& a, const int b) { return a - T(b); } 
    template<typename T> const complex<T> operator-(const int a, const complex<T>& b) { return T(a) - b; }  
    template<typename T> const complex<T> operator*(const complex<T>& a, const int b) { return a * T(b); } 
    template<typename T> const complex<T> operator*(const int a, const complex<T>& b) { return T(a) * b; }  
    template<typename T> const complex<T> operator/(const complex<T>& a, const int b) { return a / T(b); } 
    template<typename T> const complex<T> operator/(const int a, const complex<T>& b) { return T(a) / b; }  
    
    
} // namespace omega

#endif // OMEGA_FUNCTORS_HPP
