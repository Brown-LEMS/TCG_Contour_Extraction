/*
 * Live_Wire_2D_DP_gap_cpt.cpp - c++ implementation
 * This is the algorithm of Live-Wire 2D dynamic programming based on the algorithm introduced in intelligent scissors
 * This version is tailored for gap completion where boundary strengh is weak, and gradient direction not reliable
 *
 * Usage:
 *  [conf_mat, back_p_mat, len_mat] = Live_Wire_2D_DP(E_map, E_map_soft, O_map, start_point, h,w);
 *  input: E_map is [h X w] binary map, indicating exisiting contour edges
 *         E_map_soft is [h X w] boundary probability map
 *         O_map is orientation map of edges, ranging in [0, PI]
 *         start_point is [x, y, theta] of starting point
 *         h is the height
 *         w is the width
 *  output: conf_mat is [h X w] conf map
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
    double conf;
    Node(int x, int y, double theta, double c): x(x), y(y), theta(theta), conf(c){}
};

// define a comparision based on node's conf, put highest on top by default
class mycomparison
{
  bool reverse;
  public:
    mycomparison(const bool& revparam=false) {reverse=revparam;}
    
  bool operator() (const Node* n1, const Node* n2) const
  {
    if (reverse) return (n1->conf > n2->conf);
    else return (n1->conf < n2->conf);
  }
};

// this conf is directional from p to q, based on definition in 
double step_conf(Node* p, Node* q, double *E_map, double *E_map_soft, double *O_map, double *len_map, int h, int w, double* s_p)
{
    //double wZ = 0; // weight for distance transform
    double wG = 0.65; // weight for gradient magnitude
    double wS = 0.3; //  weight for shape term
    double wD = 1-wG-wS;
    int px = p->x;
    int py = p->y;
    int qx = q->x;
    int qy = q->y;
    
    /*********** fG: Reward from image gradient *******************/        
    double fG = E_map_soft[qx*h + qy];
//     fG = std::min(fG,0.3)/0.3;
    fG = fG*fG/(fG*fG+0.01);
    
    /*********** fZ: penalize the neighboorhood of existing contours *******************/        
    double fZ = E_map[qx*h + qy];
    if(fZ==1) //  do not penalize ending on the existing contour
        fZ==0;


    /*********** fO: agreement with gradient direction *******************/    
    double dx = (qx - px);
    double dy = (qy - py);
    double dist = sqrt(dx*dx + dy*dy);
    
    double dq = cos(O_map[qx*h + qy]) * dx + sin(O_map[qx*h + qy]) * dy;
    dq = dq/dist;
    if(dq <0)
        dq = -dq;

    double fO = exp(-acos(dq)*acos(dq)/M_PI*4/M_PI*4/2);
        
    /*********** fS: shape term for gap completion *******************/
    double fS = 1/M_PI * acos( (cos(p->theta)*dx + sin(p->theta)*dy )/dist);
    if(px==s_p[0] && py==s_p[1])
        wS *= 2;
    
    /*********** step reward *******************/
    // NOTE: *dist is important, or the different sampling sparcity in vertical horizontal and diagonal direction will affect  
    double step_s = 0;
    if(len_map[px*h + py]>3)
        step_s += wG*fG*(1-fZ)*fO*dist - wS*fS*dist - wD*dist;
    else
        step_s += 2*wG*fG*fO*dist - wS*fS*dist - wD*dist;
    
    return step_s;
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
    double theta_th = (mxGetPr(prhs[6]))[0];
    
    
    // create the output matrix: conf of optimal path ending at each p
    plhs[0] = mxCreateDoubleMatrix(h, w, mxREAL);
    // create a pointer to the output matrix
    double *conf_mat = mxGetPr(plhs[0]);
    
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
    my_pri_queue L(mycomparison(false)); // priority queue, keep nodes sorting in conf, max on top
    std::vector<int> ind_L ( h*w ,0); // indicator for checking if a node is in L
    std::vector<int> ind_e ( h*w ,0); // indicator for checking if a node is expended
    std::vector<Node*> All_Nodes ( h*w ,NULL); ; // full list of node holding key and conf
    for (int x=0; x<w; x++)
    {
        for (int y=0; y<h; y++)
        {
            Node* cur_n = new Node(x, y, O_map[x*h + y], 0);
            All_Nodes[x*h + y] = cur_n;
            conf_mat[x*h + y] = 0; // for safe, initialize the output
            back_p_mat[x*h + y] = 0;
            len_mat[x*h + y]=0;
        }
    }
    conf_mat[int(s[0])*h + int(s[1])] = 0;
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
                if(cos_dir_diff<0.1) //  skip the angle beyond 90 degree
                    continue;
                double dx0 = rx-s[0];
                double dy0 = ry-s[1];
                double cos_dir_diff0 = (cos(s[2])*dx0 + sin(s[2])*dy0)/sqrt(dx0*dx0 + dy0*dy0);
                if(cos_dir_diff0 < cos(theta_th)) //  skip the angle beyond search angel from starting point
                    continue;
                //std::cout << "rx:" << rx << " " << "ry:" << ry << std::endl;
                Node* r_node =  All_Nodes[rx*h+ry];
                double conf_step = step_conf(q_node, r_node, E_map, E_map_soft, O_map, len_mat, h, w, s);
                if(conf_step<=0)
                    continue;
                
                double conf_temp = q_node->conf + conf_step;
                //std::cout << "conf_temp" << conf_temp << std::endl;
                double len_temp = len_mat[qx*h+qy] + sqrt(dx*dx + dy*dy);
                //std::cout << "len_temp" <<len_temp << std::endl;
                                
                // if r_node is in L, and conf_temp > conf of r_node 
                if(ind_L[rx*h+ry] && conf_temp > r_node->conf){
                    // re-asign conf and back pointer
                    conf_mat[rx*h+ry] = conf_temp;
                    r_node->conf = conf_temp;
                    r_node->theta = atan2(dy, dx);
                    back_p_mat[rx*h+ry] = qx*h+qy;
                    len_mat[rx*h+ry] = len_temp;
                }
                // r_node not in L
                else if (!ind_L[rx*h+ry]){ 
                    // assign conf and back pointer of r_node
                    conf_mat[rx*h+ry] = conf_temp;
                    r_node->conf = conf_temp;
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