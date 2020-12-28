#include<stdio.h>
#include <math.h>
#include "mex.h"


#define INF 1e10
// 	getLU(A, b, c, d, epsilon, bestBranch, bits_per_dim, dimension, num_thread, pN, L_Bound, U_Bound);

void getLU(double *A, double *b, double *c, double *d, double epsilon, 
                   double *bestbranch, double bits_per_dim_, double dimension_, 
                   double num_thread_, int numPoint, double* L_Bound, double* U_Bound) {
    
    double getmaxorminAx( double value, double paramax, double paramin, int flag );
    double getmaxorminSquare( double valmax, double valmin, int flag );
    
    int bits_per_dim = (int)bits_per_dim_;
    int dimension = (int)dimension_;
    int num_thread = (int)num_thread_;
    
	for ( int index = 0; index < num_thread; index++ ) {
        
        int ind[dimension];
        for ( int i = 0; i < dimension; ++i ) {
            ind[i] = index / (int)(pow(bits_per_dim, i));
        }       

        double subbranch[dimension * 2];
        for ( int i = 0;  i < dimension; ++i ) {
            subbranch[i] = bestbranch[i] + (bestbranch[i+dimension] / (int)bits_per_dim) * (ind[i] - (ind[i] / (int)bits_per_dim) * (int)bits_per_dim);
            subbranch[i+dimension] = subbranch[i] + bestbranch[i+dimension] / (int)bits_per_dim;
        }

        //upper and lower bounds estimation
        double e1_min = 0.0, e1_max = 0.0, e2_min = 0.0, e2_max = 0.0, e3_min = 0.0, e3_max = 0.0;
        double r_min = 0.0;
        double rou = 0.0, rou0 = 0.0;
        double e01 = 0.0, e02 = 0.0, e03 = 0.0, e0 = 0.0;
        double l = 0, u = 0;
    
        for ( int i = 0; i < numPoint; ++i ) {
            //
            e1_max = 0.0; e1_min = 0.0; e2_max = 0.0; e2_min = 0.0; e3_max = 0.0; e3_min = 0.0;
            for ( int j = 0; j < dimension; ++j ) {
                e1_max += getmaxorminAx( A[2*dimension*i + 2*j], subbranch[j + dimension], subbranch[j], 1);
                e1_min += getmaxorminAx( A[2*dimension*i + 2*j], subbranch[j + dimension], subbranch[j], -1);
                e2_max += getmaxorminAx( A[2*dimension*i + 2*j + 1], subbranch[j + dimension], subbranch[j], 1);
                e2_min += getmaxorminAx( A[2*dimension*i + 2*j + 1], subbranch[j + dimension], subbranch[j], -1);
                e3_max += getmaxorminAx( c[dimension*i + j], subbranch[j + dimension], subbranch[j], 1);
                e3_min += getmaxorminAx( c[dimension*i + j], subbranch[j + dimension], subbranch[j], -1);                                          
            }
            e1_max += b[2*i];
            e1_min += b[2*i];
            e2_max += b[2*i + 1];
            e2_min += b[2*i + 1];            
            e3_max += d[i];            
            e3_min += d[i];    
            
           //     
            r_min = 0.0;
            r_min += getmaxorminSquare( e1_max, e1_min, -1);
            r_min += getmaxorminSquare( e2_max, e2_min, -1);
            r_min -= epsilon * epsilon * getmaxorminSquare( e3_max, e3_min, 1);
            
	        if ( r_min > 0 ) {
                rou = 1;
            }
            else {
                rou = 0;
            }
            
            l = l + rou;                
            
            e01 = 0.0; e02 = 0.0; e03 = 0.0; e0 = 0.0;
            for ( int j = 0; j < dimension; ++j ) {
                e01 += A[2*dimension*i + 2*j] * 0.5 * ( subbranch[j] + subbranch[j + dimension] );
                e02 += A[2*dimension*i + 2*j + 1] * 0.5 * ( subbranch[j] + subbranch[j + dimension] );
                e03 += c[dimension*i + j]* 0.5 * ( subbranch[j] + subbranch[j + dimension] );
            }           
            e01 += b[2*i];
            e02 += b[2*i + 1];
            e03 += d[i];
            e0 += e01 * e01 + e02 * e02 - epsilon * epsilon * e03 * e03;
            
           if (e0<0) {
                rou0 = 0;
            }
            else {
                rou0 = 1;
            }

            u = u + rou0;           
            
        }
 
		L_Bound[index] = l / numPoint;
		U_Bound[index] = u / numPoint;       
  
	}
}



double getmaxorminAx( double value, double paramax, double paramin, int flag ) {
    
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

double getmaxorminSquare( double valmax, double valmin, int flag ) {
    
    double result = 0.0;
    if ( flag == 1 ) { // max(value^2), value in [valmin, valmax]
        result = valmax * valmax >= valmin * valmin ? valmax * valmax : valmin * valmin;                
    }    
    
    if ( flag == -1 ) { // min(value^2), value in [valmin, valmax]
        if (valmax * valmin > 0) {
            result = valmax * valmax >= valmin * valmin ? valmin * valmin : valmax * valmax;                 
        }
        else {
            result = 0;
        }        
    }
    
    return result;
}




void mexFunction(int nlhs, mxArray *plhs[], 
                            int nrhs, const mxArray *prhs[])
{
    
	double *A, *b, *c, *d, *bestBranch;
	double bits_per_dim_, dimension_, num_thread_;
	double epsilon;

	A = mxGetPr(prhs[0]);
	b = mxGetPr(prhs[1]); 
	c = mxGetPr(prhs[2]); 
	d = mxGetPr(prhs[3]); 
	epsilon = *(mxGetPr(prhs[4])); 
	bestBranch = mxGetPr(prhs[5]);
	bits_per_dim_ = *(mxGetPr(prhs[6])); 
	dimension_ = *(mxGetPr(prhs[7]));
	num_thread_ = *(mxGetPr(prhs[8]));

	int pN;	
	pN = mxGetN(prhs[3]); 

	double *L_Bound, *U_Bound;

	plhs[0] = mxCreateDoubleMatrix(1, (int)num_thread_, mxREAL); //Êä³ö
	plhs[1] = mxCreateDoubleMatrix(1, (int)num_thread_, mxREAL); //Êä³ö
	L_Bound = mxGetPr(plhs[0]);
	U_Bound = mxGetPr(plhs[1]);


	getLU(A, b, c, d, epsilon, bestBranch, bits_per_dim_, dimension_, num_thread_, pN, L_Bound, U_Bound);

}








