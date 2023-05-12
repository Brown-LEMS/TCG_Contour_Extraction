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
    
	double *c3 = mxGetPr(prhs[2]);
    const int *c3_dims;
    c3_dims    = mxGetDimensions(prhs[2]);
    
    double *min_dist = mxGetPr(prhs[3]);
    
    //cout<<"min_dist: "<<min_dist[0]<<endl;
        
    int *output_dims;
    output_dims  = (int *)mxMalloc(2*sizeof(int));
    output_dims[0]=1;
    output_dims[1]=1;

	plhs[0] = mxCreateNumericArray(2, output_dims, mxDOUBLE_CLASS, mxREAL);
	double *b = mxGetPr(plhs[0]);
    int length_1 = c1_dims[1];
    double *e1 = (double *)mxMalloc(c1_dims[0]*sizeof(double));
    double *e2 = (double *)mxMalloc(c2_dims[0]*sizeof(double));
    double *e3 = (double *)mxMalloc(c3_dims[0]*sizeof(double));
    e2 = &c2[0];
    e3 = &c3[0];
    //cout<<"e2: "<<e2[0]<<" "<<e2[1]<<" "<<e2[2]<<endl;
    //cout<<"e3: "<<e3[0]<<" "<<e3[1]<<" "<<e3[2]<<endl;
    
    b[0] = 0;
    int sum_inside_edges = 0;
    
    for (int i=0; i<length_1; i++)
    {        
        e1 = &c1[i*c1_dims[0]];
        double dist_0 = compute_edgel_dist(e1, e2);
        double dist_1 = compute_edgel_dist(e1, e3);
        if(dist_0 < min_dist[0] && dist_1 < min_dist[0])
            sum_inside_edges ++;

        if (sum_inside_edges >=12 || double(sum_inside_edges)/double(length_1)>0.6)
        {
            b[0] = 1; 
            break;
        }
    }
    
    //cout<<"sum_inside_edges: "<<sum_inside_edges<<endl;

    
    mxFree(output_dims);
}