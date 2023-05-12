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
    double *cost_thresh = mxGetPr(prhs[2]);
    
    //cout<<"c1 dims "<<c1_dims[0]<<" "<<c1_dims[1]<<endl;
    //cout<<"c2 dims "<<c2_dims[0]<<" "<<c2_dims[1]<<endl;
    
    
    int *output_dims;
    output_dims  = (int *)mxMalloc(2*sizeof(int));
    output_dims[0]=1;
    output_dims[1]=1;

	plhs[0] = mxCreateNumericArray(2, output_dims, mxDOUBLE_CLASS, mxREAL);
	double *cost = mxGetPr(plhs[0]);
    int length_1 = c1_dims[1];
    int length_2 = c2_dims[1];
    double *e1  = (double *)mxMalloc(c1_dims[0]*sizeof(double));
    double *e2  = (double *)mxMalloc(c2_dims[0]*sizeof(double));
    
    for (int i=0; i<length_1; i++)
    {        
        e1 = &c1[i*c1_dims[0]];
        double min_value = 1000;
        int min_j = 0;
        for (int j=0; j<length_2; j++)
        {
            e2 = &c2[j*c2_dims[0]];
            double d = compute_edgel_dist(e1, e2);
            if(d<min_value)
            {
                min_value = d;
                min_j = j;
            }
        }
        
        e2 = &c2[min_j*c2_dims[0]];
        double d_theta = compute_ori_diff(e1, e2);
        //cout<< "e1 "<< e1[0] <<" "<<e1[1]<<" "<<e1[2]<<endl;
        //cout<< "e2 "<< e2[0] <<" "<<e2[1]<<" "<<e2[2]<<endl;
        //cout<<"min_j "<<min_j<<" min dist "<< min_value << " d_theta "<<d_theta<<endl;
        
        if(min_value < cost_thresh[0] && d_theta < (PI*0.25))
            sum_overlap ++;
        sum += (min_value + d_theta);
    }
    
    //mxFree(e1);
    //mxFree(e2);
    
	cost[0] = sum/double(length_1);
    
    if(double(sum_overlap)/double(length_1) < 0.7 && sum_overlap <20)
        cost[0] = 1000;
    mxFree(output_dims);
}


