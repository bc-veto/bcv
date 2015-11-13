#include <algorithm>
#include <iostream>
#include <string.h>
#include "mex.h"
#include "followup.hpp"
#include <lal/LALDetectors.h>

// logSkymap = wposteriors(channelNames, rate, wSw, xSw, directions, taus)

using namespace std;
using namespace omega;

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    apply f;
 
    if (nrhs != 6)
        mexErrMsgTxt("Five arguments required");
    if (nlhs != 1)
        mexErrMsgTxt("One output argument required");
    
    // channelNames
    
    if (!mxIsCell(prhs[0]))
        mexErrMsgTxt("channelNames must be a cell array");
    
	int numberOfDetectors = mxGetNumberOfElements(prhs[0]);
    
	f.detectors.resize(numberOfDetectors);
    
    for (int i = 0; i != numberOfDetectors; ++i)
    {
        mxArray* cp = mxGetCell(prhs[0], i);
        if (!mxIsChar(cp))
            mexErrMsgTxt("channelNames{i} must be char precision");
        if (mxGetNumberOfElements(cp) < 2)
            mexErrMsgTxt("channelNames{i} must have at least two characters");
        char c[3];
        c[0] = mxGetChars(cp)[0];
        c[1] = mxGetChars(cp)[1];
        c[2] = 0;
        if (strcmp(c, "T1") == 0)
        {
            f.detectors[i] = LAL_TAMA_300_DETECTOR;
        }
        else if (strcmp(c, "V1") == 0)
        {
            f.detectors[i] = LAL_VIRGO_DETECTOR;
        }
        else if (strcmp(c, "G1") == 0)
        {
            f.detectors[i] = LAL_GEO_600_DETECTOR;
        }
        else if (strcmp(c, "H2") == 0)
        {
            f.detectors[i] = LAL_LHO_2K_DETECTOR;
        }
        else if (strcmp(c, "H1") == 0)
        {
            f.detectors[i] = LAL_LHO_4K_DETECTOR;
        }
        else if (strcmp(c, "L1") == 0)
        {
            f.detectors[i] = LAL_LLO_4K_DETECTOR;
        }
        else
        {
            mexErrMsgTxt("channelNames{i} has unrecognized prefix");
        }
    }
    
    // rate
    
    if (!mxIsDouble(prhs[1]))
        mexErrMsgTxt("rate must be a double array");
    if (mxGetNumberOfElements(prhs[1]) != 1)
        mexErrMsgTxt("rate must be a scalar");
    
    f.rate = (int) mxGetScalar(prhs[1]);
    
    // wSw
    
    if (!mxIsDouble(prhs[2]))
        mexErrMsgTxt("wSw must be a double array");
    if (mxGetNumberOfElements(prhs[2]) != numberOfDetectors)
        mexErrMsgTxt("wSw must have same number of elements as channelNames");
 	f.wSw.resize(numberOfDetectors);
    copy(mxGetPr(prhs[2]), mxGetPr(prhs[2]) + numberOfDetectors, f.wSw.begin());
    
    // xSw
 
	if (!mxIsCell(prhs[3]))
		mexErrMsgTxt("xSw must be a cell array");
	if (mxGetNumberOfElements(prhs[3]) != numberOfDetectors)
        mexErrMsgTxt("xSw must have same number of cells as channelNames");
    f.xSw_real.resize(numberOfDetectors);
    f.xSw_imag.resize(numberOfDetectors);
    for (int i = 0; i != numberOfDetectors; ++i)
    {
        mxArray* cp = mxGetCell(prhs[3], i);
        if (!mxIsDouble(cp))
            mexErrMsgTxt("xSw{i} must be a double array");
        double* p = mxGetPr(cp);
        int n = mxGetNumberOfElements(cp);
        f.xSw_real[i].resize(n);
        copy(p, p + n, f.xSw_real[i].begin());
        p = mxGetPi(cp);
        if (!p)
            mexErrMsgTxt("xSw{i} must have an imaginary components");
        f.xSw_imag[i].resize(n);
        copy(p, p + n, f.xSw_imag[i].begin());
    }
    
    // directions
    
    if (!mxIsDouble(prhs[4]))
        mexErrMsgTxt("directions must be a double array");
    if (mxGetM(prhs[4]) != 2)
        mexErrMsgTxt("directions must be a 2 x N matrix");
    int n = mxGetN(prhs[4]);
    f.directions.resize(n);
    for (int i = 0; i != n; ++i)
    {
        double* p = mxGetPr(prhs[4]);
        f.directions[i].x[0] = p[i * 2];
        f.directions[i].x[1] = p[i * 2 + 1];
    }
    
    // taus
    
    if (!mxIsDouble(prhs[5]))
        mexErrMsgTxt("taus must be a double array");
    f.taus.resize(mxGetNumberOfElements(prhs[5]));
    copy(mxGetPr(prhs[5]), mxGetPr(prhs[5]) + f.taus.size(), f.taus.begin());
    
	// perform the followup
    
    f();
    
    plhs[0] = mxCreateDoubleMatrix(1, f.directions.size(), mxREAL);
    double* p = mxGetPr(plhs[0]);
    
    copy(f.log_skymap.begin(), f.log_skymap.end(), p);
    
}
