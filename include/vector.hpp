#ifndef OMEGA_VECTOR_HPP
#define OMEGA_VECTOR_HPP

#include <algorithm>

#include "functors.hpp"

namespace omega
{
    
    using namespace std;
     
    namespace detail
    {
        
        class heap_vector {};
        
        template<size_t N> class stack_vector {};
        
        template <class T, class A> class range_vector {};
        template <class T, class A> class const_range_vector {};     
        
        template<class T, class A, class F> class unary_mapped {};       
        
        class index_vector {};
        class zero_vector {};
        template<class A, class B> class join_vector {};
        template<class A> class zero_pad_vector {};        
        
    } // namespace detail

    template<typename T, typename A = detail::heap_vector> class vector;
    
    // computed assignment operators
    
    template<class T, class A, class U> vector<T, A>& operator+=(vector<T, A>& a, const U& b)
    {
        for (size_t i = 0; i != a.size(); ++i)
            a[i] += b;
        return a;
    }
    
    template<class T, class A, class U, class B> vector<T, A>& operator+=(vector<T, A>& a, const vector<U, B>& b)
    {
        for (size_t i = 0; i != a.size(); ++i)
            a[i] += b[i];
        return a;
    }   
    
    template<class T, class A, class U> vector<T, A>& operator-=(vector<T, A>& a, const U& b)
    {
        for (size_t i = 0; i != a.size(); ++i)
            a[i] -= b;
        return a;
    }
    
    template<class T, class A, class U, class B> vector<T, A>& operator-=(vector<T, A>& a, const vector<U, B>& b)
    {
        for (size_t i = 0; i != a.size(); ++i)
            a[i] -= b[i];
        return a;
    }     
    
    template<class T, class A, class U> vector<T, A>& operator*=(vector<T, A>& a, const U& b)
    {
        for (size_t i = 0; i != a.size(); ++i)
            a[i] *= b;
        return a;
    }
    
    template<class T, class A, class U, class B> vector<T, A>& operator*=(vector<T, A>& a, const vector<U, B>& b)
    {
        for (size_t i = 0; i != a.size(); ++i)
            a[i] *= b[i];
        return a;
    }     
    
    template<class T, class A, class U> vector<T, A>& operator/=(vector<T, A>& a, const U& b)
    {
        for (size_t i = 0; i != a.size(); ++i)
            a[i] /= b;
        return a;
    }
    
    template<class T, class A, class U, class B> vector<T, A>& operator/=(vector<T, A>& a, const vector<U, B>& b)
    {
        for (size_t i = 0; i != a.size(); ++i)
            a[i] /= b[i];
        return a;
    }          
    
    // non-elementwise functions
    
    template<class T, class A> const T sum(const vector<T, A>& a)
    {
        T b = T();
        for (size_t i = 0; i != a.size(); ++i)
            b += a[i];
        return b;
    }
    
    template<class T, class A> const T min(const vector<T, A>& a)
    {
        assert(a.size() > 0);
        T b = a[0];
        for (size_t i = 1; i != a.size(); ++i)
            if (a[i] < b)
                b = a[i];
        return b;
    }
    
    template<class T, class A> const size_t argmax(const vector<T, A>& a)
    {
        assert(a.size() > 0);
        size_t b = 0;
        for (size_t i = 1; i != a.size(); ++i)
            if (a[b] < a[i])
                b = i;
        return b;
    }
    
    template<class T, class A> const T max(const vector<T, A>& a)
    {
        assert(a.size() > 0);
        T b = a[0];
        for (size_t i = 1; i != a.size(); ++i)
            if (b < a[i])
                b = a[i];
        return b;
    }      
    
    template<class T, class A> T mean(const vector<T, A>& a)
    {
        return sum(a) / a.size();
    }   
    
    template<class T, class A> const vector<T, A> cumsum(const vector<T, A>& a)
    {
        assert(a.size() > 0);
        vector<T> b(a.size());
        b[0] = a[0];
        for (size_t i = 1; i != a.size(); ++i)
            b[i] = b[i - 1] + a[i];
        return b;
    }        
    
    // subscripting functionality
    
    class range
    {
        range(); // prohibit
    public:
        size_t a_, b_;
        range(size_t a, size_t b) : a_(a), b_(b) { assert(a_ <= b_); }
    };
    
    template<typename T> class vector<T, detail::heap_vector>
    {
        T* _;
        size_t n;
        
    public:
        
        // utility function
        
        void swap(vector& b) { std::swap(_, b._); std::swap(n, b.n); }
        
        // basic semantics
        
        vector() : _(0), n(0) {}
        
        vector (const vector& a) 
            : _(a.size() ? new T[a.size()] : 0), n(a.size()) 
        { 
            for (size_t i = 0; i != n; ++i) 
                _[i] = a[i]; 
        }
        
        ~vector() { delete[] _; }       
        
        vector& operator=(const vector& a) 
        { 
            vector b(a); 
            swap(b); 
            return *this; 
        }       

        // sized construction
        
        explicit vector(const size_t i) : _(i ? new T[i] : 0), n(i) {} 
        template<typename U> vector(const size_t i, const U& a) : _(i ? new T[i] : 0), n(i) 
        {
            fill(begin(), end(), a);
        }
        
        // conversions

        template<typename U, typename B> vector(const vector<U, B>& b) 
            : _(b.size() ? new T[b.size()] : 0), n(b.size()) 
        { 
            for (size_t i = 0; i != n; ++i)
                _[i] = b[i];
        }        

        template<typename U, typename B> 
        vector& operator=(const vector<U, B>& a)
        {
            if (n == a.size())
            {
                for (size_t i = 0; i != n; ++i)
                    _[i] = a[i];
            }
            else
            {
                vector b(a);
                swap(b);
            }
            return *this;
        }
        
        // size
        
        const size_t size() const { return n; }        
        void resize(const size_t i) { vector b(i); swap(b); }
        
        // iterators
        
        T* begin() { return _; }
        T* end() { return _ + n; }
        const T* begin() const { return _; }
        const T* end() const { return _ + n; }
        
        // subscripting
        
        T& operator[](const size_t i) { assert(i < n); return _[i]; }
        const T& operator[](const size_t i) const { assert(i < n); return _[i]; }
        
        vector<T, detail::range_vector<T, detail::heap_vector> > operator[](const range& r) 
        { 
            assert(r.b_ <= size());
            return vector<T, detail::range_vector<T, detail::heap_vector> >(r, *this);
        }
        
        vector<T, detail::const_range_vector<T, detail::heap_vector> > operator[](const range& r) const
        { 
            assert(r.b_ <= size());
            return vector<T, detail::const_range_vector<T, detail::heap_vector> >(r, *this);
        }
                
    }; // heap vector
    
    template<class T, class A> vector<T> evaluate(const vector<T, A>& a)
    {
        return a;
    }
    
    template<class T, size_t N> class vector<T, detail::stack_vector<N> >
    {
        T _[N];
    public:
        template<class U, class B> vector(const vector<U, B>& b)
        {
            assert(b.size() == N);
            for (size_t i = 0; i != N; ++i)
                _[i] = b[i];            
        }
        template<class U, class B> vector& operator=(const vector<U, B>& b)
        {
            assert(b.size() == N);
            for (size_t i = 0; i != N; ++i)
                _[i] = b[i];
        }
        const size_t size() const { return N; }
        T* begin() { return _; }
        T* end() { return _ + N; }
        const T* begin() const { return _; }
        const T* end() const { return _ + N; }       
        
        // subscripting
        
        T& operator[](const size_t i) { assert(i < N); return _[i]; }
        const T& operator[](const size_t i) const { assert(i < N); return _[i]; }
        
        vector<T, detail::range_vector<T, detail::stack_vector<N> > > operator[](const range& r) 
        { 
            assert(r.b_ <= size());
            return vector<T, detail::range_vector<T, detail::stack_vector<N> > >(r, *this);
        }
        
        vector<T, detail::const_range_vector<T, detail::stack_vector<N> > > operator[](const range& r)  const
        { 
            assert(r.b_ <= size());
            return vector<T, detail::const_range_vector<T, detail::stack_vector<N> > >(r, *this);
        }
        
    }; // stack vector
        
    // expression templates for the form f(vector)
          
    template<class T, class A> class vector<T, detail::range_vector<T, A > >
    {
        const range& r_;
        vector<T, A>& a_;
        vector();
    public:
        vector(const range& r, vector<T, A>& a) : r_(r), a_(a) {}    
        vector& operator=(const vector& a)
        {
            assert((r_.b_ - r_.a_) == a.size());
            for (size_t i = 0; i != a.size(); ++i)
                a_[r_.a_ + i] = a[i];
            return *this;
        }
        template<class U, class B> vector& operator=(const vector<U, B>& a)
        {
            assert((r_.b_ - r_.a_) == a.size());
            for (size_t i = 0; i != a.size(); ++i)
                a_[r_.a_ + i] = a[i];
            return *this;
        }
        
        const size_t size() const { return r_.b_ - r_.a_; }
        const T operator[](const size_t i) const
        {
            assert(r_.a_ + i < r_.b_);
            return a_[r_.a_ + i];
        }
        T& operator[](const size_t i)
        {
            assert(r_.a_ + i < r_.b_);
            return a_[r_.a_ + i];
        }    
        
    }; // range vector
    
    template<class T, class A> class vector<T, detail::const_range_vector<T, A> >
    {
        const range& r_;
        const vector<T, A>& a_;
        vector();
    public:
        vector(const range& r, const vector<T, A>& a) : r_(r), a_(a) {}
        const size_t size() const { return r_.b_ - r_.a_; }
        const T operator[](const size_t i) const
        {
            assert(r_.a_ + i < r_.b_);
            return a_[r_.a_ + i];
        }
    }; // const range vector
    
    // functions mapping reals to reals, or complex to complex

    template<class T, class U, class A, class F> class vector<T, detail::unary_mapped<U, A, F> >
    {
        const vector<U, A>& _;
        const F f_;
    public:
        vector(const vector<U, A>& a, const F& f) : _(a), f_(f) {}
        const T operator[](size_t i) const 
        {
            assert(i < _.size());
            return f_(_[i]); 
        }
        size_t size() const { return _.size(); }
    };     
    
    #define UNARY_VECTOR_OPERATION(NAME) \
    template<class T, class A> \
        const vector<typename NAME##_functor<T>::type, detail::unary_mapped<T, A, NAME##_functor<T> > > \
        NAME(const vector<T, A>& a) \
    { \
        return vector<typename NAME##_functor<T>::type, detail::unary_mapped<T, A, NAME##_functor<T> > >(a, NAME##_functor<T>()); \
    }    
    
    UNARY_VECTOR_OPERATION(abs)
    UNARY_VECTOR_OPERATION(acos)
    UNARY_VECTOR_OPERATION(asin)
    UNARY_VECTOR_OPERATION(atan)
    UNARY_VECTOR_OPERATION(cos)
    UNARY_VECTOR_OPERATION(conj)
    UNARY_VECTOR_OPERATION(cosh)
    UNARY_VECTOR_OPERATION(exp)
    UNARY_VECTOR_OPERATION(log)
    UNARY_VECTOR_OPERATION(log10)
    UNARY_VECTOR_OPERATION(sin)
    UNARY_VECTOR_OPERATION(sinh)
    UNARY_VECTOR_OPERATION(sqrt)
    UNARY_VECTOR_OPERATION(tan)
    UNARY_VECTOR_OPERATION(tanh) 
    UNARY_VECTOR_OPERATION(real) 
    UNARY_VECTOR_OPERATION(imag) 
    UNARY_VECTOR_OPERATION(arg) 
    UNARY_VECTOR_OPERATION(norm) 
    
    #undef UNARY_VECTOR_OPERATION
        
    // expression templates for the form f(vector, vector)
        
    template<class T, class A, class U, class B, template<class, class> class F> class binary_mapped {};
    
    template<class S, class T, class A, class U, class B, template <class, class> class F> class vector<S, binary_mapped<T, A, U, B, F> >
    {
        const vector<T, A>& a_;
        const vector<U, B>& b_;
    public:
        explicit vector(const vector<T, A>& a, const vector<U, B>& b) : a_(a), b_(b) 
        {
            assert(a.size() == b.size());
        }
        const S operator[](size_t i) const 
        {
            assert(i < a_.size());
            return F<T, U>()(a_[i], b_[i]); 
        }
        size_t size() const { return a_.size(); }
    };
    
    #define BINARY_VECTOR_OPERATION(NAME, OPERATOR)\
     template<class T, class A, class U, class B> \
        vector<typename NAME##_functor<T, U>::type, binary_mapped<T, A, U, B, NAME##_functor> > \
            operator OPERATOR (const vector<T, A>& a, const vector<U, B>& b) \
    { \
        return vector<typename NAME##_functor<T, U>::type, binary_mapped<T, A, U, B, NAME##_functor> >(a, b); \
    } \
    \
    template<class T, class A, class U> \
        vector<typename bind2nd<NAME##_functor<T, U>, T, U>::type, detail::unary_mapped<T, A, bind2nd<NAME##_functor<T, U>, T, U> > > \
            operator OPERATOR (const vector<T, A>& t, const U& u) \
    { \
        return \
            vector<typename bind2nd<NAME##_functor<T, U>, T, U>::type, detail::unary_mapped<T, A, bind2nd<NAME##_functor<T, U>, T, U> > > \
                (t, bind2nd<NAME##_functor<T, U>, T, U>(NAME##_functor<T, U>(), u)); \
    } \
    \
    template<class T, class U, class B> \
        vector<typename bind1st<NAME##_functor<T, U>, T, U>::type, detail::unary_mapped<U, B, bind1st<NAME##_functor<T, U>, T, U> > > \
            operator OPERATOR (const T& t, const vector<U, B>& u) \
    { \
        return \
            vector<typename bind1st<NAME##_functor<T, U>, T, U>::type, detail::unary_mapped<U, B, bind1st<NAME##_functor<T, U>, T, U> > > \
                (u, bind1st<NAME##_functor<T, U>, T, U>(NAME##_functor<T, U>(), t)); \
    }
    
    BINARY_VECTOR_OPERATION(plus, +)
    BINARY_VECTOR_OPERATION(minus, -)
    BINARY_VECTOR_OPERATION(times, *)
    BINARY_VECTOR_OPERATION(divided_by, /)
    
    #undef BINARY_VECTOR_OPERATION
    
    #define BINARY_VECTOR_FUNCTION(NAME)\
     template<class T, class A, class U, class B> \
        vector<typename NAME##_functor<T, U>::type, binary_mapped<T, A, U, B, NAME##_functor> > \
            NAME(const vector<T, A>& a, const vector<U, B>& b) \
    { \
        return vector<typename NAME##_functor<T, U>::type, binary_mapped<T, A, U, B, NAME##_functor> >(a, b); \
    } \
    \
    template<class T, class A, class U> \
        vector<typename bind2nd<NAME##_functor<T, U>, T, U>::type, detail::unary_mapped<T, A, bind2nd<NAME##_functor<T, U>, T, U> > > \
            NAME(const vector<T, A>& t, const U& u) \
    { \
        return \
            vector<typename bind2nd<NAME##_functor<T, U>, T, U>::type, detail::unary_mapped<T, A, bind2nd<NAME##_functor<T, U>, T, U> > > \
                (t, bind2nd<NAME##_functor<T, U>, T, U>(NAME##_functor<T, U>(), u)); \
    } \
    \
    template<class T, class U, class B> \
        vector<typename bind1st<NAME##_functor<T, U>, T, U>::type, detail::unary_mapped<U, B, bind1st<NAME##_functor<T, U>, T, U> > > \
            NAME(const T& t, const vector<U, B>& u) \
    { \
        return \
            vector<typename bind1st<NAME##_functor<T, U>, T, U>::type, detail::unary_mapped<U, B, bind1st<NAME##_functor<T, U>, T, U> > > \
                (u, bind1st<NAME##_functor<T, U>, T, U>(NAME##_functor<T, U>(), t)); \
    }  
    
    BINARY_VECTOR_FUNCTION(pow)
    BINARY_VECTOR_FUNCTION(atan2)
    
    #undef BINARY_VECTOR_FUNCTION      
            
    template<> class vector<size_t, detail::index_vector> 
    {
        size_t n_;
        vector();
        vector& operator=(const vector&);
    public:
        explicit vector(const size_t n) : n_(n) {}
        const size_t operator[](const size_t i) const
        {
            assert(i < size());
            return i;
        }
        const size_t size() const { return n_; }
    };
    
    typedef vector<size_t, detail::index_vector> indices;
        
    template<class T> class vector<T, detail::zero_vector>
    {
        size_t n_;
        vector();
    public:
        explicit vector(const size_t n) : n_(n) {}
        const size_t size() const { return n_; }
        const T operator[](const size_t i) const
        {
            assert(i < size());
            return 0;
        }
    };
    
    template<typename T> const vector<T, detail::zero_vector> zeros(const size_t n)
    {
        return vector<T, detail::zero_vector>(n);
    }
           
    template<class T, class A, class B> class vector<T, detail::join_vector<A, B> >
    {
        const vector<T, A>& a_;
        const vector<T, B>& b_;
    public:
        vector(const vector<T, A>& a, const vector<T, B>& b) : a_(a), b_(b) {}
        const size_t size() const { return a_.size() + b_.size(); }
        const T operator[](const size_t i) const
        {
            assert(i < size());
            return i < a_.size() ? a_[i] : b_[i - a_.size()];
        }
    };
    
    template<class T, class A, class B> const vector<T, detail::join_vector<A, B> > join(const vector<T, A>& a, const vector<T, B>& b)
    {
        return vector<T, detail::join_vector<A, B> >(a, b);
    }
    

    template< class T, class A> class vector<T, detail::zero_pad_vector<A> >
    {
        const size_t a_;
        const vector<T, A>& b_;
        const size_t c_;
    public:
        vector(const size_t a, const vector<T, A>& b, const size_t c) : a_(a), b_(b), c_(c) 
        {
            assert(a_ >= 0);
            assert(c_ >= 0);
        }
        const size_t size() const { return a_ + b_.size() + c_; }
        const T operator[](const size_t i) const
        {
            assert(i < size());
            return i < a_ ? 0 : (i < a_ + b_.size() ? b_[i - a_] : 0);
        }
        
    };
    
    template<class T, class A> 
    const vector<T, detail::zero_pad_vector<A> > 
        zero_pad(const size_t a, const vector<T, A>& b, const size_t c)
    {
        return vector<T, detail::zero_pad_vector<A> >(a, b, c);
    }
        
    // not elementwise, but require temporary vector
       
    template<class T, class A> const T median(const vector<T, A>& a)
    {
        vector<T> b(a);
        nth_element(&b[0], &b[b.size() / 2], &b[b.size()]);
        return b[b.size() / 2];
    }
    
    namespace detail
    {
        template<class T, class A> class rancor
        {
            const vector<T, A>& a_;
        public:
            explicit rancor(const vector<T, A>& a) : a_(a) {}
            bool operator()(const size_t i, const size_t j) const { return a_[i] > a_[j]; }
        };   
    }
    
    // index of greatest, next greatest, ... least elements
    template<class T, class A> const vector<size_t> rank(const vector<T, A>& a)
    {

        vector<size_t> b(indices(a.size()));
        std::sort(b.begin(), b.end(), detail::rancor<T, A>(a));
        return b;
    }      
    
}

#endif // OMEGA_VECTOR_HPP
