
#ifndef __match_theta_hh__
#define __match_theta_hh__

class Matrix;

// returns the cost of the assignment
double matchEdgeMaps (
    const Matrix& bmap1, const Matrix& bmap2,
    const Matrix& theta1, const Matrix& theta2,
    double theta_diff,double maxDist, double outlierCost,
    Matrix& match1, Matrix& match2);

#endif // __match_theta_hh__
