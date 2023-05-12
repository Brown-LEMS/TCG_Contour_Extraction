#include "util_functions.h"
#include "mex.h"
using namespace std;

void mexFunction( int nlhs, mxArray *plhs[] , int nrhs, const mxArray *prhs[] )
{
    
    int sum_overlap = 0;
    double sum = 0;
    
    
	double *c1 = mxGetPr(prhs[0]);
    const int *c1_dims;
    c1_dims    = mxGetDimensions(prhs[0]);
	double *c2 = mxGetPr(prhs[1]);
    const int *c2_dims;
    c2_dims    = mxGetDimensions(prhs[1]);
    double *local_dist = mxGetPr(prhs[2]);
    
    //cout<<"c1 dims "<<c1_dims[0]<<" "<<c1_dims[1]<<endl;
    //cout<<"c2 dims "<<c2_dims[0]<<" "<<c2_dims[1]<<endl;
    
    
    int *output_dims;
    output_dims  = (int *)mxMalloc(2*sizeof(int));
    output_dims[0]=1;
    output_dims[1]=1;

	plhs[0] = mxCreateNumericArray(2, output_dims, mxDOUBLE_CLASS, mxREAL);
	double *b = mxGetPr(plhs[0]);
    int length_1 = c1_dims[1];
    int length_2 = c2_dims[1];
    double *e1  = (double *)mxMalloc(c1_dims[0]*sizeof(double));
    double *e2  = (double *)mxMalloc(c2_dims[0]*sizeof(double));
    b[0] = 0;
    
    if (length_2 >=5) // only consider curve with more than five edges
    {
        for (int i=0; i<length_2; i++)
        {        
            e2 = &c2[i*c2_dims[0]];
            for (int j=0; j<length_1; j++)
            {
                e1 = &c1[j*c1_dims[0]];
                double d = compute_edgel_dist(e1, e2);
                if(d<2*local_dist[0] && compute_ori_diff(e1,e2) < 0.16*PI)
                {
                    sum_overlap ++;
                    break;
                }
            }

            if (sum_overlap >=12 || double(sum_overlap)/double(length_2)>0.6)
            {
                b[0] = 1; 
                break;
            }
        }
    }
    
    mxFree(output_dims);
}


