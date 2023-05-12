#include <math.h>
#include <iostream>
#define PI M_PI

double compute_edgel_dist(double *e1, double *e2)
{
	return sqrt((e1[0]-e2[0])*(e1[0]-e2[0])+(e1[1]-e2[1])*(e1[1]-e2[1]));
}

double compute_ori_diff(double *e1, double *e2)
{
	double dir_1 = e1[2];
	double dir_2 = e2[2];

	if(dir_1 >= PI)
		dir_1 -= PI;
	if(dir_2 >= PI)
		dir_2 -= PI;

	double d_theta = fabs(dir_2-dir_1);
	if(d_theta >= PI*0.5)
		d_theta = PI - d_theta;
	return d_theta;
}