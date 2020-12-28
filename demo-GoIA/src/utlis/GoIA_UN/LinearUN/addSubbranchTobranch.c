#include "mex.h"
#include "math.h"
#include<stdio.h>

void addTobranch(double *L, double *U, double best_U, double *bestBranch, 
                            double bits_per_dim, double dimension, double *branch, double numthread, double* num) {

    double* findChild(double *bestbranch, int best_U_index, int bits_per_dim,  int dimension);
	double* branchM ;
	branchM = (double *)mxMalloc(sizeof(double)*(int)(2 * dimension + 2));
    for(int i=0; i<numthread; i++){
		//printf("%lf\n", L[i]);

        if(L[i] <= best_U){

            if( (L[i]>-10)&&(U[i]<10) ){
                branchM = findChild(bestBranch, i, (int)bits_per_dim, (int)dimension);
                branchM[2*(int)dimension] = L[i];
                branchM[2*(int)dimension + 1] = U[i];
                
				for (int j = 0; j < (int)(2 * dimension + 2); j++){
					branch[j + (int)num[0]*(int)(2 * dimension + 2)] = branchM[j];
				}					

                num[0] = num[0]+1;
				//printf("%lf\n", num[0]);
            }
        }
        
   }
}



    double* findChild(double *bestbranch, int best_U_index, int bits_per_dim, int dimension)
    {
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
        
        double* best_U_branch;
        best_U_branch = (double *)mxMalloc(sizeof(double)*(int)(2 * dimension + 2));
        
        for (int i = 0; i < dimension; ++i) {
            best_U_branch[i] = subbranch[i];
            best_U_branch[i+dimension] = width[i];
        }

        return best_U_branch;

    }




void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{

    double *L, *U, *bestBranch;
    double best_U;
    double bits_per_dim;
    double dimension;
    double numthread;

    L = mxGetPr(prhs[0]);
    U = mxGetPr(prhs[1]);
    best_U = *mxGetPr(prhs[2]);
    bestBranch = mxGetPr(prhs[3]);
    bits_per_dim = *mxGetPr(prhs[4]);
    dimension = *mxGetPr(prhs[5]);  
    numthread = *mxGetPr(prhs[6]);

    double* num;
    plhs[0] = mxCreateDoubleMatrix((int)(2*dimension+2), (int)numthread, mxREAL);
    double* branch;
    branch = mxGetPr(plhs[0]);
    plhs[1] = mxCreateDoubleMatrix(1, 1, mxREAL);
    num = mxGetPr(plhs[1]);
    num[0] = 0;

    addTobranch(L, U, best_U, bestBranch, bits_per_dim, dimension, branch, numthread, num);

//printf("%lf", num);
//
//
//printf("%lf\n", branch[2]);
//mxFree(branch);

}
