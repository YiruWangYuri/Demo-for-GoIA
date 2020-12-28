#include<stdio.h>
#include <math.h>
#include "mex.h"


#define INF 1e10
// 	getLU(data, m_data, n_data, epsilon, bestBranch, m_bB, n_bB, bits_per_dim, dimension, num_thread, L_Bound, U_Bound);

void getLU(double *data, int numPoint, double epsilon, 
                  double *bestbranch, double bits_per_dim_, double dimension_, 
                  double num_thread_,  double *L_Bound, double *U_Bound) {

    double getmaxormin(double value, double paramax, double paramin, int flag);
    int bits_per_dim = (int)bits_per_dim_;
    int dimension = (int)dimension_;
    int num_thread = (int)num_thread_;
    
    for( int index = 0; index < num_thread; index++) {
        
        int ind[dimension];
        for ( int i = 0; i < dimension; ++i ) {
            ind[i] = index / (int)(pow(bits_per_dim, i));
        }

        double subbranch[dimension * 2];
        for ( int i = 0;  i < dimension; ++i) {
            subbranch[i] = bestbranch[i] + (bestbranch[i+dimension] / (int)bits_per_dim) * (ind[i] - (ind[i] / (int)bits_per_dim) * (int)bits_per_dim);
            subbranch[i+dimension] = subbranch[i] + bestbranch[i+dimension] / (int)bits_per_dim;
        }
        
   
        //selectBranch
        double r2max = 0.0, r2min = 0.0;
        
        for ( int i = 0; i < dimension; ++i ) {
            
            if ( subbranch[i + dimension] * subbranch[i + dimension]  > subbranch[i] * subbranch[i] )                                                                           {
                r2max += ( subbranch[i + dimension] * subbranch[i + dimension] );
                r2min += ( (subbranch[i] * subbranch[i + dimension]) >0 ) * subbranch[i] * subbranch[i];
            }
            else {
                r2max += ( subbranch[i] * subbranch[i] );
                r2min += ( (subbranch[i] * subbranch[i + dimension]) >0 ) * subbranch[i + dimension] * subbranch[i + dimension];
            }
            
        }
            
        if (r2max < 1)	{
            L_Bound[index] = -INF;
            U_Bound[index] = +INF;
            continue;
        }         
        
        if (r2min > 1)	{
            L_Bound[index] = -INF;
            U_Bound[index] = +INF;
            continue;
        }           

    //upper and lower bounds estimation
        
        double e_max = 0.0, e_min = 0.0, e0 = 0.0;
        double r_min = 0.0, r_max = 0.0;
        double rou = 0.0, rou0 = 0.0;
        double l = 0.0, u = 0.0;
        double tempbranchCenter[dimension];
        double tempCenterNorm = 0.0;
        
        for ( int i = 0; i < numPoint; i++ ) {            
            e_max = 0.0; e_min = 0.0;
            for ( int j = 0; j < dimension; ++j ) {
                e_max += getmaxormin(data[i + j*numPoint], subbranch[j + dimension], subbranch[j], 1);
                e_min += getmaxormin(data[i + j*numPoint], subbranch[j + dimension], subbranch[j], -1);
            }
//             e_max += subbranch[dimension*2-1];
//             e_min += subbranch[dimension - 1];


            if ( e_max * e_min<0 ) {
                r_min = 0;
                e_max = fabs(e_max);
                e_min = fabs(e_min);
                //r_max = max(e_max,e_min);
                if ( e_max>e_min ) {
                    r_max = e_max;
                }
                else {
                    r_max = e_min;
                }
            }
            else {
                e_max = fabs(e_max);
                e_min = fabs(e_min);
                //r_min = min(e_max,e_min);
                if ( e_max>e_min ) {
                    r_min = e_min;
                }
                else {
                    r_min = e_max;
                }
                //r_max = max(e_max,e_min);
                if ( e_max>e_min ) {
                    r_max = e_max;
                }
                else {
                    r_max = e_min;
                }
            }

            if ( r_min > epsilon ) {
                rou = 1;
            }
            else {
                rou = 0;
            }
            
            l = l + rou;     
        

            tempCenterNorm = 0.0;
            for ( int j = 0; j < dimension; ++j ) { 
                tempbranchCenter[j] = 0.5 * (subbranch[j] + subbranch[j + dimension]);
                tempCenterNorm += tempbranchCenter[j] * tempbranchCenter[j];
            }
            tempCenterNorm = sqrt(tempCenterNorm);
            for ( int j = 0; j < dimension; ++j ) { 
                tempbranchCenter[j] /= tempCenterNorm;
            }
            
            e0 = 0.0;
            for ( int j = 0; j < dimension; ++j ) {
                e0 +=  tempbranchCenter[j] * data[i + j*numPoint];
            }
            e0 = fabs(e0);
 
            if ( e0 < epsilon ) {
                rou0 = 0;
            }
            else {
                rou0 = 1;
            }

            u = u + rou0;
	}



    L_Bound[index] = l/numPoint;
    U_Bound[index] = u/numPoint;

   }

}


void mexFunction(int nlhs, mxArray *plhs[], 
                            int nrhs, const mxArray *prhs[]) {

    size_t m_data, n_data, m_bB, n_bB;
	double *data, *bestBranch;
	double bits_per_dim_, dimension_, num_thread_;
	double epsilon;
    
    m_data = mxGetM(prhs[0]); // get the row number of input data
    n_data = mxGetN(prhs[0]); // get the column number of input data
    m_bB = mxGetM(prhs[2]); // get the row number of input bestBrach
    n_bB = mxGetN(prhs[2]); // get the row number of input data
    
    data = mxGetPr(prhs[0]);   
	epsilon = *(mxGetPr(prhs[1]));
	bestBranch = mxGetPr(prhs[2]); 
	bits_per_dim_ = *(mxGetPr(prhs[3]));
	dimension_ = *(mxGetPr(prhs[4]));
	num_thread_ = *(mxGetPr(prhs[5]));

	double *L_Bound, *U_Bound;

	plhs[0] = mxCreateDoubleMatrix(1, (int)num_thread_, mxREAL); 
	plhs[1] = mxCreateDoubleMatrix(1, (int)num_thread_, mxREAL); 
	L_Bound = mxGetPr(plhs[0]);
	U_Bound = mxGetPr(plhs[1]);

	getLU(data, m_data, epsilon, bestBranch, bits_per_dim_, dimension_, num_thread_, L_Bound, U_Bound);

}


double getmaxormin(double value, double paramax, double paramin, int flag) {
    
    double result = 0.0;
    
    if ( flag == 1 ) { // max(para*value)

        if ( value >= 0 )
            { result = paramax*value; }
        else
            { result = paramin*value; }
                
    }
    
    if (flag == -1) {// min(para*value)
        
        if (value >= 0)
            { result = paramin*value; }
        else
            { result = paramax*value; }
        
    }
    
    return result;

}





