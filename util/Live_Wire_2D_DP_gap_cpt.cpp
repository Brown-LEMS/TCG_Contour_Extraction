/*
 * Live_Wire_2D_DP_gap_cpt.cpp - c++ implementation
 * This is the algorithm of Live-Wire 2D dynamic programming based on the algorithm introduced in intelligent scissors
 * This version is tailored for gap completion where boundary strengh is weak, and gradient direction not reliable
 *
 * Usage:
 *  [cost_mat, back_p_mat, len_mat] = Live_Wire_2D_DP(E_map, E_map_soft, O_map, start_point, h,w);
 *  input: E_map is [h X w] binary map, indicating exisiting contour edges
 *         E_map_soft is [h X w] boundary probability map
 *         O_map is orientation map of edges, ranging in [0, PI]
 *         start_point is [x, y, theta] of starting point
 *         h is the height
 *         w is the width
 *  output: cost_mat is [h X w] cost map
 *          back_p_mat is [h X w] map recording the back pointer which is coded in x*h + y
 *          len_mat is [h X w] map recording the length of optimal path
 *
 * This is a MEX file for MATLAB.
 */

#include <stdlib.h>
#include <math.h>
#include "matrix.h"
#include "mex.h"
#include <iostream>
#include <queue>          // std::priority_queue
#include <vector>         // std::vector
#include <functional>     // std::greater

struct Node {
    int x;
    int y;
    double theta;
    double cost;
    Node(int x, int y, double theta, double c): x(x), y(y), theta(theta), cost(c){}
};

// define a comparision based on node's cost, put highest on top by default
class mycomparison
{
  bool reverse;
  public:
    mycomparison(const bool& revparam=false) {reverse=revparam;}
    
  bool operator() (const Node* n1, const Node* n2) const
  {
    if (reverse) return (n1->cost > n2->cost);
    else return (n1->cost < n2->cost);
  }
};

// this cost is directional from p to q, based on definition in 
double step_cost(Node* p, Node* q, double *E_map, double *E_map_soft, double *O_map, int h, int w, double* s_p)
{
    double wZ = 0.2; // weight for zero crossing
    double wG = 0.6; // weight for gradient magnitude
    double wD = 0; // weight for gradient agreement
    double wS = 1-wZ-wG-wD; //  weight for shape term
    int px = p->x;
    int py = p->y;
    int qx = q->x;
    int qy = q->y;
    double fZ = E_map[qx*h + qy];
    double fG = 1-E_map_soft[qx*h + qy];
    fG = std::min(fG,0.7);

    double dx = (qx - px);
    double dy = (qy - py);
    double dist = sqrt(dx*dx + dy*dy);
    
    
    double dp;
    // starting point is in [-pi, pi], need to be enforced
    if(px==s_p[0] && py==s_p[1])
        dp = cos(s_p[2]) * dx + sin(s_p[2]) * dy;
    else    
        dp = cos(O_map[px*h + py]) * dx + sin(O_map[px*h + py]) * dy;
    dp = dp/dist;
    double dq = cos(O_map[qx*h + qy]) * dx + sin(O_map[qx*h + qy]) * dy;
    dq = dq/dist;

    // such that dp is always positive
    // orientation map consider bi-directions
    // starting point is in [-pi, pi], need to be enforced
    if(dp <0 && !(px==s_p[0] && py==s_p[1] ))
    {
        dp = -dp;
    }
    if(dq <0)
        dq = -dq;

    //double fD = 1/M_PI * (acos(dp) + acos(dq));
    double fD = exp(-acos(dq)*acos(dq)/M_PI*6/M_PI*6);
        
    /*********** fS shape term for gap completion *******************/
    double fS = 1/M_PI * acos( (cos(p->theta)*dx + sin(p->theta)*dy )/dist);
    if(px==s_p[0] && py==s_p[1])
        wS *= 5;
    
    // NOTE: *dist is important, or the different sampling sparcity in vertical horizontal and diagonal direction will affect    
    //return wZ*fZ*dist + wD*fD*dist + wG*fG*dist + wS*fS*dist;
    //return wZ*fZ*dist + wG*fG*fD*dist + wS*fS*dist;
    return wZ*fZ*dist + wG*(1- (1-fG)*fD)*dist + wS*fS*dist;
}

/*  the gateway routine.  */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {

    /****************** variables declaration/initialization ******************/

    //  edges map binary
    double *E_map = mxGetPr(prhs[0]);
    // edge map soft
    double *E_map_soft = mxGetPr(prhs[1]);
    // orienatation map
    double *O_map = mxGetPr(prhs[2]);
    // start point: (x, y, theta)
    double *s = mxGetPr(prhs[3]);
    // size of the input maps 
    int h = (mxGetPr(prhs[4]))[0];
    int w = (mxGetPr(prhs[5]))[0];
    
    
    // create the output matrix: cost of optimal path ending at each p
    plhs[0] = mxCreateDoubleMatrix(h, w, mxREAL);
    // create a pointer to the output matrix
    double *cost_mat = mxGetPr(plhs[0]);
    
    // create the output matrix: save back pointer 
    plhs[1] = mxCreateDoubleMatrix(h, w, mxREAL);
    // create a pointer to the output matrix
    double *back_p_mat = mxGetPr(plhs[1]);

    // create the output matrix: save length of optimal path 
    plhs[2] = mxCreateDoubleMatrix(h, w, mxREAL);
    // create a pointer to the output matrix
    double *len_mat = mxGetPr(plhs[2]);
    
    //std::cout << "fetch input and create output" << std:: endl;
    
    /****************** variables declaration/initialization ******************/
    // using mycomparison:
    typedef std::priority_queue<Node*, std::vector<Node*>, mycomparison> my_pri_queue;
    my_pri_queue L(mycomparison(true)); // priority queue, keep nodes sorting in cost, withi min on top
    std::vector<int> ind_L ( h*w ,0); // indicator for checking if a node is in L
    std::vector<int> ind_e ( h*w ,0); // indicator for checking if a node is expended
    std::vector<Node*> All_Nodes ( h*w ,NULL); ; // full list of node holding key and cost
    for (int x=0; x<w; x++)
    {
        for (int y=0; y<h; y++)
        {
            Node* cur_n = new Node(x, y, O_map[x*h + y], 0);
            All_Nodes[x*h + y] = cur_n;
            cost_mat[x*h + y] = DBL_MAX; // for safe, initialize the output matrix
            back_p_mat[x*h + y] = 0;
            len_mat[x*h + y]=0;
        }
    }
    cost_mat[int(s[0])*h + int(s[1])] = 0;
    len_mat[int(s[0])*h + int(s[1])] = 0;
    All_Nodes[int(s[0])*h + int(s[1])]->theta = s[2];
    // it might cause seg fault to change input pointer, the enforce is done in temp_compute anyway
    // O_map[int(s[0])*h + int(s[1])] =  s[2]; //  enforce the orientation of starting point

    //std::cout << "complete initialization of variables" << std:: endl;

    
    /****************** dynamic programming start ******************/
    L.push(All_Nodes[s[0]*h + s[1]]);
    ind_L[int(s[0])*h + int(s[1])] = 1;
    while(!L.empty()){
        Node* q_node = L.top();
        L.pop();
        int qx = q_node->x;
        int qy = q_node->y;
        double qtheta = q_node->theta;
        //std::cout << "qx:" << qx << " qy:" << qy << " w:" << w <<" h:" << h << " sx:" << s[0] << " sy:" << s[1] << std::endl;
        ind_L[qx*h+qy] = 0;
        ind_e[qx*h+qy] = 1;
        
        // for all the neighboring nodes, that not expended
        for (int dx = -1; dx <=1; dx++){
            int rx =  qx + dx;
            if(rx<0 || rx >=w)
                continue;
            for(int dy= -1; dy<=1; dy++){
                int ry =  qy + dy;
                if(ry<0 || ry >=h)
                    continue;
                if(rx == qx && ry == qy) // skip the center
                    continue;
                if(ind_e[rx*h+ry]) // skip the expended node
                    continue;
                double cos_dir_diff = (cos(qtheta)*dx + sin(qtheta)*dy)/sqrt(dx*dx + dy*dy);
                if(cos_dir_diff<0) //  skip the angle beyond 90 degree
                    continue;
                
                //std::cout << "rx:" << rx << " " << "ry:" << ry << std::endl;
                Node* r_node =  All_Nodes[rx*h+ry];
                double cost_temp = q_node->cost + step_cost(q_node, r_node, E_map, E_map_soft, O_map, h, w, s);
                //std::cout << "cost_temp" << cost_temp << std::endl;
                double len_temp = len_mat[qx*h+qy] + sqrt(dx*dx + dy*dy);
                //std::cout << "len_temp" <<len_temp << std::endl;
                
                // if r_node is in L, and cost_temp < cost of r_node 
                if(ind_L[rx*h+ry] && cost_temp < r_node->cost){
                    // re-asign cost and back pointer
                    cost_mat[rx*h+ry] = cost_temp;
                    r_node->cost = cost_temp;
                    r_node->theta = atan2(dy, dx);
                    back_p_mat[rx*h+ry] = qx*h+qy;
                    len_mat[rx*h+ry] = len_temp;
                }
                // r_node not in L
                else if (!ind_L[rx*h+ry]){ 
                    // assign cost and back pointer of r_node
                    cost_mat[rx*h+ry] = cost_temp;
                    r_node->cost = cost_temp;
                    r_node->theta = atan2(dy, dx);
                    back_p_mat[rx*h+ry] = qx*h+qy;
                    len_mat[rx*h+ry] = len_temp;
                    /************ additional constraints for gap completion: existing contour should not be any other node's parent *********/
                    // if reach a existing contour point, do not insert into into L
                    // this makes the search terminate much earlier in many cases
                    if(E_map[rx*h+ry]==1) 
                    {
                        ind_e[rx*h+ry]=1;
                        continue;
                    }
                    // push r_node into L
                    L.push(r_node);
                    ind_L[rx*h+ry] = 1;
                }

            }
        }
    }
    
    // clear memory
    for (int x=0; x<w; x++)
    {
        for (int y=0; y<h; y++)
        {
            Node* cur_n = All_Nodes[x*h + y];
            delete cur_n;
        }
    }
    //std::cout << "done" << std::endl;
    
    
}