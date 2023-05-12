#include "util_functions.h"
#include "mex.h"
#include <vector>

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
	double *cost = mxGetPr(plhs[0]);
    int length_1 = c1_dims[1];
    int length_2 = c2_dims[1];
    double *e1  = (double *)mxMalloc(c1_dims[0]*sizeof(double));
    double *e2  = (double *)mxMalloc(c2_dims[0]*sizeof(double));
    
    vector<double> v_1(length_1, 1000), v_2(length_2, 1000);
    // make two vectors to make the matched edges and unmateched edges
    for (int i=0; i<length_1; i++)
    {        
        e1 = &c1[i*c1_dims[0]];
        double min_dist = 1000;
        int min_j = 0;
        for (int j=0; j<length_2; j++)
        {
            e2 = &c2[j*c2_dims[0]];
            double d = compute_edgel_dist(e1, e2);
            if(d<min_dist)
            {
                min_dist = d;
                min_j = j;
            }
        }
        e2 = &c2[min_j*c2_dims[0]];
        if(min_dist<local_dist[0])
            v_1[i]=pow(min_dist/local_dist[0], 4) + 0.6*compute_ori_diff(e1, e2);
    }
    
    for (int i=0; i<length_2; i++)
    {        
        e2 = &c2[i*c2_dims[0]];
        double min_dist = 1000;
        int min_j = 0;
        for (int j=0; j<length_1; j++)
        {
            e1 = &c1[j*c1_dims[0]];
            double d = compute_edgel_dist(e1, e2);
            if(d<min_dist)
            {
                min_dist = d;
                min_j = j;
            }
        }
        e1 = &c1[min_j*c1_dims[0]];
        if(min_dist<local_dist[0])
            v_2[i]=pow(min_dist/local_dist[0], 4) + 0.6*compute_ori_diff(e1, e2);
    }
    
    //compute the transformation cost as the length of unmatched {1000, 1000, ...1000} fragments plus the total localization dist
	double deform_length_1 = 0, deform_length_2 = 0, local_cost_1 = 0, local_cost_2 = 0;
    for (int i = 0; i< (v_1.size()-1);i++)
    {
        e1 = &c1[i*c1_dims[0]];
        e2 = &c1[(i+1)*c1_dims[0]];
        double dist = compute_edgel_dist(e1, e2);
        if(dist > 10)
            dist = 10;
        if(v_1[i]==1000 && v_1[i+1]==1000)
            deform_length_1 += dist*local_dist[0];
        if(v_1[i]!=1000 && v_1[i+1]!=1000)
            local_cost_1 += (v_1[i]+v_1[i+1])*dist/2;
    }
    
    for (int i = 0; i< (v_2.size()-1);i++)
    {
        e1 = &c2[i*c2_dims[0]];
        e2 = &c2[(i+1)*c2_dims[0]];
        double dist = compute_edgel_dist(e1, e2);
        if(dist > 10)
            dist = 10;
        if(v_2[i]==1000 && v_2[i+1]==1000)
            deform_length_2 += dist*local_dist[0];
        if(v_2[i]!=1000 && v_2[i+1]!=1000)
            local_cost_2 += (v_2[i]+v_2[i+1])*dist/2;
    }
    
    cost[0] = deform_length_1 + deform_length_2 + (local_cost_1 + local_cost_2)/2;
    
    mxFree(output_dims);
}


