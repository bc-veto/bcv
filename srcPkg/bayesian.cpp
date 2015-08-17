// bayesian takes double-whitened strain data from the Omega pipeline,
// performs a bayesian follow-up analysis on a specified t-f region, and
// returns parameter estimates and robust detection statistics

// C++ standard library

#include <algorithm>
#include <complex>
#include <iostream>
#include <limits>
#include <list>
#include <stdexcept>
#include <set>

// FFTW

#include <fftw3.h>

// LAL

#include <lal/LALConstants.h>
#include <lal/Random.h>
#include <lal/Skymap.h>

// Omega

#include "interval.hpp"
#include "vector.hpp"
//#include "exp_double.hpp"
#include "followup.hpp"

namespace omega {
    
    using namespace std;

    // we need estimates of the variation in K and x over the parameter ranges
    // abstract these out so we can try different methodologies, i.e. interval
    // computation, point sampling etc.
    
    // also abstract out bisection (keep it to bisection)
 
    ostream& operator<<(ostream& a, const interval& x)
    {
        return a << '[' << x.a << ", " << x.b << ']';
    }    
                
    template<int N> void chunk<N>::compute(apply& a, XLALSkymapPlanType& plan)
    {
        // get ranges of elements of K
        // get ranges of times
        // get ranges of strains

        interval delay[N]; 
        interval k[N][N];
        interval logNormalization;

        for (int i = 0; i != 2; ++i)
        {
            double theta = i ? x[0].b : x[0].a;
            for (int j = 0; j != 2; ++j)
            {
                double phi = i ? x[1].b : x[1].a;
                double tp[2] = { theta, phi };
                XLALSkymapDirectionPropertiesType properties;
                XLALSkymapDirectionPropertiesConstruct(&plan, tp, &properties);
                XLALSkymapKernelType kernel;
                XLALSkymapKernelConstruct(&plan, &properties, a.wSw.begin(), &kernel);
                if ((i == 0) && (j == 0))
                {
                    for (int l = 0; l != N; ++l)
                    {
                        delay[l] = properties.delay[l];
                        for (int m = 0; m != N; ++m)
                        {
                            k[l][m] = kernel.k[l][m];
                        }
                    }
                    logNormalization = kernel.logNormalization;
                }
                else
                {
                    for (int l = 0; l != N; ++l)
                    {
                        delay[l] = hull(delay[l], properties.delay[l]);
                        for (int m = 0; m != N; ++m)
                        {
                            k[l][m] = hull(k[l][m], kernel.k[l][m]);
                        }
                    }
                    logNormalization = hull(logNormalization, kernel.logNormalization);
                }
            }
        }

        //for (int i = 0; i != N; ++i)
        //{
         //   for (int j = 0; j != N; ++j)
         //   {
         //       cout << k[i][j];
        //    }            
        //    cout << endl;
        //}
        
        interval s[N];
        interval xr[N];
        interval xi[N];
        
        //cout << length(x[2]) * a.rate << endl;
        //cout << x[2].a << " " << x[2].b << endl;
        
        for (int i = 0; i != N; ++i)
        {
            s[i] = (x[2] + delay[i]) * interval(plan.sampleFrequency);
            int b = floor(s[i].a);
            int e = ceil(s[i].b) + 1;
            
            //cout << x[i] << endl;
                        
            xr[i] = interval(
                    *min_element(a.xSw_real[i].begin() + b, a.xSw_real[i].begin() + e),
                    *max_element(a.xSw_real[i].begin() + b, a.xSw_real[i].begin() + e));
            xi[i] = interval(
                    *min_element(a.xSw_imag[i].begin() + b, a.xSw_imag[i].begin() + e),
                    *max_element(a.xSw_imag[i].begin() + b, a.xSw_imag[i].begin() + e));
            
            
            
            cout << e - b << xr[i] << xi[i] << endl;
            //cout << e - b << endl;
        }
        
        logp = interval(0.);
        for (int i = 0; i != N; ++i)
        {
            for (int j = 0; j != N; ++j)
            {
                logp += (xr[i] * xr[j] + xi[i] * xi[j]) * k[i][j];
            }
        }
        logp = logp * interval(.5) + logNormalization * interval(2.);
        
        logf = logp + interval(log(length(x[0]) * length(x[1]) * length(x[2]))) + log(sin(x[0]));
        
        logError = log1p(-exp(logf.a - logf.b)) + logf.b;
        
    }
    
    template<int N> bool operator<(const chunk<N>& a, const chunk<N>& b)
    {
        return a.logError > b.logError;
    }
    
    template<int N> void apply::adaptive_n()
    {
        XLALSkymapPlanType plan;
        XLALSkymapPlanConstruct(rate, N, detectors.begin(), &plan);
        
        double min_t;
        double max_t;
        
        min_t = taus[detectors.size() * 2];
        max_t = taus[detectors.size() * 2 + 1];
        
        //cout << min_t << ", " << max_t << endl;

        multiset<chunk<N> > working;
        
        chunk<N> c;
        c.x[0] = interval(0, pi);
        c.x[1] = interval(0, 2 * pi);
        c.x[2] = interval(min_t, max_t);
        c.compute(*this, plan);
        working.insert(c);
        
        bool flag = true;
        
        while (flag)
        {
            chunk<N> c = *working.begin();
            working.erase(working.begin());

            //cout << length(c.x[0]) * 180/pi << endl;
            //cout << c.logp << endl;
            
            for (int i = 0; i != 2; ++i)
            {
                chunk<N> d;
                d.x[0] = partition(c.x[0], i, 2);
                for (int j = 0; j != 2; ++j)
                {
                    d.x[1] = partition(c.x[1], j, 2);
                    for (int k = 0; k != 2; ++k)
                    {
                        d.x[2] = partition(c.x[2], k, 2);
                        d.compute(*this, plan);
                        working.insert(d);
                    }
                }
            }
            
            //interval t;
            //for (typename multiset<chunk<N> >::iterator i = working.begin(); i != working.end(); ++i)
            //{
            //    t += exp(i->logf);
            //}
            //cout << t << endl;
            
            //flag = working.size() < 10;
        }
        
                 
        
    };
    
    void apply::adaptive()
    {
        switch (detectors.size())
        {
        case 1:
            adaptive_n<1>();
            break;
        case 2:
            adaptive_n<2>();
            break;
        case 3:
            adaptive_n<3>();
            break;
        }
    }
    
   
    void test()
    {     
        
        RandomParams* rng;
        rng = XLALCreateRandomParams(0);
        
        int siteNumbers[] = { LAL_LHO_4K_DETECTOR, LAL_LLO_4K_DETECTOR, LAL_VIRGO_DETECTOR, LAL_GEO_600_DETECTOR, LAL_LHO_2K_DETECTOR };                
        XLALSkymapPlanType plan;
        XLALSkymapPlanConstruct(16384, 3, siteNumbers, &plan);
        
        double direction[2];
        direction[0] = LAL_PI * XLALUniformDeviate(rng);
        direction[1] = LAL_TWOPI * XLALUniformDeviate(rng);

        XLALSkymapDirectionPropertiesType properties;
        XLALSkymapDirectionPropertiesConstruct(&plan, direction, &properties);
                
        double wSw[5] = { 10., 10., 10., 10., 10. };
        XLALSkymapKernelType kernel;
        XLALSkymapKernelConstruct(&plan, &properties, wSw, &kernel);
        
        //cout << "Kernel:" << endl;
        //for (int i = 0; i != 3; ++i)
        //{
        //    for (int j = 0; j != 3; ++j)
        //    {
        //        cout << kernel.k[i][j] << ' ';
        //    }
        //    cout << endl;
        //}
        
        double* xSw_real[3];
        double* xSw_imag[3];

        for (int i = 0; i != 3; ++i)
        {
            xSw_real[i] = new double[plan.sampleFrequency];
            xSw_imag[i] = new double[plan.sampleFrequency];
            for (int j = 0; j != plan.sampleFrequency; ++j)
            {
                double ar = XLALNormalDeviate(rng) * sqrt(wSw[i]);
                double ai = XLALNormalDeviate(rng) * sqrt(wSw[i]);
                xSw_real[i][j] = ar;
                xSw_imag[i][j] = ai;
            }
            xSw_real[i][plan.sampleFrequency / 2] += XLALNormalDeviate(rng) * wSw[i];
            xSw_imag[i][plan.sampleFrequency / 2] += XLALNormalDeviate(rng) * wSw[i];
            //cout << ar << ' ' << ai << endl;
        }

        double realLogPosterior;
        XLALSkymapApply(&plan, &properties, &kernel, xSw_real, 0.5, &realLogPosterior);
        double imagLogPosterior;
        XLALSkymapApply(&plan, &properties, &kernel, xSw_imag, 0.5, &imagLogPosterior);
        
        cout << "Reference = " << (realLogPosterior + imagLogPosterior) << endl;

        ///////////////////////////////////////////////////////////////////
        
        apply a;
        a.detectors.resize(3);
        copy(siteNumbers, siteNumbers + 3, a.detectors.begin());
        a.rate = plan.sampleFrequency;
        a.wSw.resize(3);
        copy(wSw, wSw + 3, a.wSw.begin());
        a.xSw_real.resize(3);
        a.xSw_imag.resize(3);
        for (int i = 0; i != 3; ++i)
        {
            a.xSw_real[i].resize(a.rate);
            a.xSw_imag[i].resize(a.rate);
            copy(xSw_real[i], xSw_real[i] + a.rate, a.xSw_real[i].begin());
            copy(xSw_imag[i], xSw_imag[i] + a.rate, a.xSw_imag[i].begin());
            //for (int j = 0; j != 8; ++j)
            //{
            //    cout << xSw_real[i][j] << " ? " << a.xSw_real[i][j] << endl;
            //    cout << xSw_imag[i][j] << " ? " << a.xSw_imag[i][j] << endl;
            //}
        }
        
        chunk<3> c;
        c.x[0] = interval(direction[0]);
        c.x[1] = interval(direction[1]);
        c.x[2] = interval(0.5);
        c.compute(a, plan);
        
        cout << "Computed = " << c.logp << endl;
                
        ///////////////////////////////////////////////////////////////////
        
        a.taus.resize(3 * 2 + 3);
        a.taus[3*2] = 0.485;
        a.taus[3*2 + 1] = 0.515;
        
        a.adaptive();
        
        for (int i = 0; i != 3; ++i)
        {
             delete xSw_real[i];
             delete xSw_imag[i];
             
        }
        
    }
    
    
    // old implementation /////////////////////////////////////////////////
    
    void apply::operator()() {
        XLALSkymapPlanType plan;
        
        // construct the network
        XLALSkymapPlanConstruct(rate, detectors.size(), detectors.begin(), &plan);
        
        // allocate space for the output
        log_skymap.resize(directions.size());
        
        // unpack the timing information
        double min_t;
        double max_t;
        double delta_t;
        
        vector<double> min_ts(detectors.size());
        vector<double> max_ts(detectors.size());
        for (int i = 0; i != detectors.size(); ++i) {
            min_ts[i] = taus[i * 2];
            max_ts[i] = taus[i * 2 + 1];
        }
        min_t = taus[detectors.size() * 2];
        max_t = taus[detectors.size() * 2 + 1];
        delta_t = taus[detectors.size() * 2 + 2];
        
        int n = static_cast<int>(ceil((max_t - min_t) / delta_t));
        
        vector<double> p_t(n);
        
        vector<double*> xSw_real_p(xSw_real.size());
        vector<double*> xSw_imag_p(xSw_imag.size());
        
        for (int i = 0; i != xSw_real.size(); ++i) {
            xSw_real_p[i] = xSw_real[i].begin();
            xSw_imag_p[i] = xSw_imag[i].begin();
        }
        
        for (int i = 0; i != directions.size(); ++i) {
            XLALSkymapDirectionPropertiesType properties;
            double thetaphi[2] = { directions[i].x[0], directions[i].x[1] };
            XLALSkymapDirectionPropertiesConstruct(&plan, thetaphi, &properties);
            
            double t_begin = min_t;
            double t_end = max_t;
            
            for (int j = 0; j != detectors.size(); ++j) {
                t_begin = std::max(t_begin, min_ts[j] - properties.delay[j]);
                t_end = std::min(t_end, max_ts[j] - properties.delay[j]);
            }
            
            if (t_begin < t_end) {
                
                XLALSkymapKernelType kernel;
                XLALSkymapKernelConstruct(&plan, &properties, wSw.begin(), & kernel);
                int j = 0;
                
                for (double t = t_begin; t < t_end; t += delta_t) {
                    double log_p_real, log_p_imag;
                    XLALSkymapApply(&plan, &properties, &kernel, xSw_real_p.begin(), t, &log_p_real);
                    XLALSkymapApply(&plan, &properties, &kernel, xSw_imag_p.begin(), t, &log_p_imag);
                    p_t[j] = log_p_real + log_p_imag;
                    ++j;
                }
                
                log_skymap[i] = XLALSkymapLogTotalExp(p_t.begin(), p_t.begin() + j) - log(n);
                
            }
            else {
                // the times of interest in each detector exclude this
                // direction
                log_skymap[i] = log(0);
            }
        }
        
    }
        
    
} //namespace omega

int main(int argc, char** argv)
{
    omega::test();
    return 0;
}