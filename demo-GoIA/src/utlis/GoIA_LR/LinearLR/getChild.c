#include "mex.h"
#include "math.h"
#include<stdio.h>
void findChild(double *bestbranch, double best_U_index, double bits_per_dim,  double *best_U_branch, 
                        double dimension_) {
    
    int dimension = (int)dimension_;
    int ind[dimension];
    for ( int i = 0; i < dimension; ++i ) {
        ind[i] = best_U_index / (int)(pow((float)bits_per_dim, i));
    }

    double subbranch[dimension * 2], width[dimension];
    for ( int i = 0;  i < dimension; ++i) {
        subbranch[i] = bestbranch[i] + (bestbranch[i+dimension] / (int)bits_per_dim) * (ind[i] - (ind[i] / (int)bits_per_dim) * (int)bits_per_dim);
        width[i] = bestbranch[i+dimension] / (int)bits_per_dim;
        subbranch[i+dimension] = subbranch[i] + width[i];
    }
        
    for (int i = 0; i < dimension; ++i) {
        best_U_branch[i] = subbranch[i];
        best_U_branch[i+dimension] = width[i];
    }
   
}




void mexFunction(int nlhs,mxArray *plhs[], int nrhs,const mxArray *prhs[]) {
//[best_U_branch]=getChild(bestBranch, best_U_index, bits_per_dim, dimension);

    double *bestBranch;
    double best_U_index;
    double bits_per_dim;
    double dimension;

    double *best_U_branch;

    bestBranch = mxGetPr(prhs[0]);
    best_U_index = *mxGetPr(prhs[1]);
    bits_per_dim = *mxGetPr(prhs[2]);
    dimension = *mxGetPr(prhs[3]);

    plhs[0] = mxCreateDoubleMatrix((int)(2*dimension), 1, mxREAL);

    best_U_branch = mxGetPr(plhs[0]);

    best_U_index = best_U_index-1;
    findChild(bestBranch, best_U_index, bits_per_dim, best_U_branch, dimension);

}
