
#include <string.h>
#include <mex.h>
#include "Matrix.hh"
#include "csa.hh"
#include "match_theta.hh"

extern "C" {

static const double maxDistDefault = 0.0075;
static const double outlierCostDefault = 100;
static const double maxAngleDistDefault =  0.523598775598299;

void 
mexFunction (
    int nlhs, mxArray* plhs[],
    int nrhs, const mxArray* prhs[])
{
    // check number of arguments
    if (nlhs < 2) {
        mexErrMsgTxt("Too few output arguments.");
    }
    if (nlhs > 4) {
        mexErrMsgTxt("Too many output arguments.");
    }
    if (nrhs < 4) {
        mexErrMsgTxt("Too few input arguments.");
    }
    if (nrhs > 7) {
        mexErrMsgTxt("Too many input arguments.");
    }

    // get arguments
    double* bmap1 = mxGetPr(prhs[0]);
    double* bmap2 = mxGetPr(prhs[1]);
    double* theta1 = mxGetPr(prhs[2]);
    double* theta2 = mxGetPr(prhs[3]);

    const double maxAngleDist = 
        (nrhs>4) ? mxGetScalar(prhs[4]) : maxAngleDistDefault;
    const double maxDist = 
        (nrhs>5) ? mxGetScalar(prhs[5]) : maxDistDefault;
    const double outlierCost = 
        (nrhs>6) ? mxGetScalar(prhs[6]) : outlierCostDefault;

    // check arguments
    if (mxGetM(prhs[0]) != mxGetM(prhs[1]) 
        || mxGetN(prhs[0]) != mxGetN(prhs[1])) {
        mexErrMsgTxt("bmap1 and bmap2 must be the same size");
    }

    // check arguments
    if (mxGetM(prhs[2]) != mxGetM(prhs[3]) 
        || mxGetN(prhs[2]) != mxGetN(prhs[3])) {
        mexErrMsgTxt("bmap1 and bmap2 must be the same size");
    }

    if (maxAngleDist < 0) {
        mexErrMsgTxt("maxAngleDist must be >= 0");
    }

    if (maxDist < 0) {
        mexErrMsgTxt("maxDist must be >= 0");
    }
    if (outlierCost <= 1) {
        mexErrMsgTxt("outlierCost must be >1");
    }
    //  mexPrintf("Angle dist: %f\n",maxAngleDist);

    // do the computation
    const int rows = mxGetM(prhs[0]);
    const int cols = mxGetN(prhs[0]);
    const double idiag = sqrt( rows*rows + cols*cols );
//    mexPrintf(" dist: %f\n",maxDist*idiag);
    const double oc = outlierCost*maxDist*idiag;
    Matrix m1, m2;
    const double cost = matchEdgeMaps(
        Matrix(rows,cols,bmap1), Matrix(rows,cols,bmap2),
        Matrix(rows,cols,theta1),Matrix(rows,cols,theta2),
        maxAngleDist,
        maxDist*idiag, oc,
        m1, m2);
    
    // set output arguments
    plhs[0] = mxCreateDoubleMatrix(rows, cols, mxREAL);
    plhs[1] = mxCreateDoubleMatrix(rows, cols, mxREAL);
    double* match1 = mxGetPr(plhs[0]);
    double* match2 = mxGetPr(plhs[1]);
    memcpy(match1,m1.data(),m1.numel()*sizeof(double));
    memcpy(match2,m2.data(),m2.numel()*sizeof(double));
    if (nlhs > 2) { plhs[2] = mxCreateScalarDouble(cost); }
    if (nlhs > 3) { plhs[3] = mxCreateScalarDouble(oc); }
}

}; // extern "C"
