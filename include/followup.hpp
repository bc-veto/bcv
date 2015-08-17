#ifndef FOLLOWUP_HPP
#define FOLLOWUP_HPP

#include <list>

#include "interval.hpp"
#include "vector.hpp"

#include <lal/Skymap.h>

namespace omega
{

	struct thetaphi
	{
		double x[2];
	};
	
	class apply
	{
        
    public:
        
        vector<int> detectors;                 // LAL_DETECTOR codes
        
        int rate;
        
        vector<double> wSw;               // Filter normalizations
        vector<vector<double> > xSw_real; // Time-domain filtered data
        vector<vector<double> > xSw_imag; // Time-domain filtered data
        
        vector<thetaphi> directions; // Directions to sample at
        vector<double> taus; // Times to sample at
        
        vector<double> log_skymap;
        
        void operator()();
        
        void adaptive();
        template<int N> void adaptive_n();
	};
    
    template<int N> struct chunk
    {
        interval x[3]; // region
        interval logp; // range of values
        interval logf; // contribution to integral
        double logError; // contribution to error
        void compute(apply& a, XLALSkymapPlanType& plan);
    };
    
} // namespace omega

#endif // FOLLOWUP_HPP
