#include "util_functions.h"
#include "mex.h"

using namespace std;

void mexFunction( int nlhs, mxArray *plhs[] , int nrhs, const mxArray *prhs[] )
{
    double sum = 0;
    double *c1 = mxGetPr(prhs[0]);
    const int *c1_dims;
    c1_dims    = mxGetDimensions(prhs[0]);
    
    int *output_dims;
    output_dims  = (int *)mxMalloc(2*sizeof(int));
    output_dims[0]=1;
    output_dims[1]=1;
    
    plhs[0] = mxCreateNumericArray(2, output_dims, mxDOUBLE_CLASS, mxREAL);
	double *cost = mxGetPr(plhs[0]);
    int length_1 = c1_dims[1];
    
    double *e1  = (double *)mxMalloc(c1_dims[0]*sizeof(double));
    double *e2  = (double *)mxMalloc(c1_dims[0]*sizeof(double));
    
    for (int i=0; i<length_1-1; i++)
    {        
        e1 = &c1[i*c1_dims[0]];
        int min_j = 0;
        e2 = &c1[(i+1)*c1_dims[0]];
        double d = compute_edgel_dist(e1, e2);
        sum += d;
    }
    
    cost[0] = sum;
    mxFree(output_dims);
}