#include<stdio.h>
#include <math.h>
#include "mex.h"


#define INF 1e10
	//getLU(x,y,z,sigma,bestbranch,bits_per_dim,dimension,num_thread,num_thread_per_block,L_Bound, U_Bound);

void getLU(double * x, double * y, double * z, double sigma, double * bestbranch, double bits_per_dim, double dimension, double num_thread, int numPoint,  double* L_Bound, double* U_Bound){
double getmaxormin(double x, double paramax, double paramin, int flag);	
    for(int index = 0;index<num_thread;index++){

		int ind1 = index / (int)(pow(bits_per_dim, 0));
		int ind2 = index / (int)(pow(bits_per_dim, 1));
		int ind3 = index / (int)(pow(bits_per_dim, 2));
		int ind4 = index / (int)(pow(bits_per_dim, 3));

	double a_min = bestbranch[0] + (bestbranch[4] / (int)bits_per_dim)*(ind1 - (ind1 / (int)bits_per_dim)*(int)bits_per_dim);
	double b_min = bestbranch[1] + (bestbranch[5] / (int)bits_per_dim)*(ind2 - (ind2 / (int)bits_per_dim)*(int)bits_per_dim);
	double c_min = bestbranch[2] + (bestbranch[6] / (int)bits_per_dim)*(ind3 - (ind3 / (int)bits_per_dim)*(int)bits_per_dim);
	double d_min = bestbranch[3] + (bestbranch[7] / (int)bits_per_dim)*(ind4 - (ind4 / (int)bits_per_dim)*(int)bits_per_dim);
	double wa = bestbranch[4] / (int)bits_per_dim;
	double wb = bestbranch[5] / (int)bits_per_dim;
	double wc = bestbranch[6] / (int)bits_per_dim;
	double wd = bestbranch[7] / (int)bits_per_dim;
	double a_max = a_min + wa;
	double b_max = b_min + wb;
	double c_max = c_min + wc;
	double d_max = d_min + wd;
    

	//selectBranch
	double a2max, a2min;
	if (a_max*a_max>a_min*a_min)	{
		a2max = a_max*a_max;
		a2min = a_min*a_min;
	}
	else	{
		a2max = a_min*a_min;
		a2min = a_max*a_max;
	}

	double b2max, b2min;
	if (b_max*b_max>b_min*b_min)	{
		b2max = b_max*b_max;
		b2min = b_min*b_min;
	}
	else	{
		b2max = b_min*b_min;
		b2min = b_max*b_max;
	}

	double c2max, c2min;
	if (c_max*c_max>c_min*c_min)	{
		c2max = c_max*c_max;
		c2min = c_min*c_min;
	}
	else	{
		c2max = c_min*c_min;
		c2min = c_max*c_max;
	}

	double d2max, d2min;
	if (d_max*d_max>d_min*d_min)	{
		d2max = d_max*d_max;
		d2min = d_min*d_min;
	}
	else	{
		d2max = d_min*d_min;
		d2min = d_max*d_max;
	}


	double r2max = a2max + b2max + c2max + d2max;

	if (r2max<1)	{
		L_Bound[index] = -INF;
		U_Bound[index] = +INF;
		continue;
	}

	double r2min = ((a_max*a_min)>0)*a2min + ((b_max*b_min)>0)*b2min + ((c_max*c_min)>0)*c2min + ((d_max*d_min)>0)*d2min;

	if (r2min>1)	{
		L_Bound[index] = -INF;
		U_Bound[index] = +INF;
		continue;
	}

	//估计上下界

	int i = 0; //索引点
	double e_max, e_min, e0;
	double r_min, r_max;
	double rou = 0, rou0 = 0;
	double l = 0, u = 0;
    double tempa, tempb, tempc, tempd, tempa_un, tempb_un, tempc_un, tempd_un;


	for (i = 0; i<numPoint; i++)
	{
        double temp1, temp2, temp3;
        temp1 = getmaxormin(x[i],a_max,a_min,1);
        temp2 = getmaxormin(y[i],b_max,b_min,1);
        temp3 = getmaxormin(z[i],c_max,c_min,1);
		e_max = temp1 + temp2 + temp3 + d_max;
        
        temp1 = getmaxormin(x[i],a_max,a_min,-1);
        temp2 = getmaxormin(y[i],b_max,b_min,-1);
        temp3 = getmaxormin(z[i],c_max,c_min,-1);        
		e_min = temp1 + temp2 + temp3 + d_min;

		if (e_max*e_min<0)		{
			r_min = 0;
			e_max = fabs(e_max);
			e_min = fabs(e_min);
			//r_max = max(e_max,e_min);
			if (e_max>e_min)	{
				r_max = e_max;
			}
			else				{
				r_max = e_min;
			}
		}
		else		{
			e_max = fabs(e_max);
			e_min = fabs(e_min);
			//r_min = min(e_max,e_min);
			if (e_max>e_min)		{
				r_min = e_min;
			}
			else			{
				r_min = e_max;
			}

			//r_max = max(e_max,e_min);
			if (e_max>e_min)			{
				r_max = e_max;
			}
			else			{
				r_max = e_min;
			}
		}



		if (r_min>sigma)		{
			rou = 1;
		}
		else{
			rou = 0;
		}
		l = l + rou;

        // u
        tempa = 0.5*(a_min+a_max); tempb = 0.5*(b_min+b_max); tempc = 0.5*(c_min+c_max); tempd = 0.5*(d_min+d_max); 
        tempa_un = tempa/sqrt(tempa*tempa+tempb*tempb+tempc*tempc+tempd*tempd);
        tempb_un = tempb/sqrt(tempa*tempa+tempb*tempb+tempc*tempc+tempd*tempd);
        tempc_un = tempc/sqrt(tempa*tempa+tempb*tempb+tempc*tempc+tempd*tempd);
        tempd_un = tempd/sqrt(tempa*tempa+tempb*tempb+tempc*tempc+tempd*tempd);
        
		e0 = fabs(tempa_un*x[i] + tempb_un*y[i] + tempc_un*z[i] + tempd_un);

		if (e0<sigma)		{
			rou0 = 0;
		}
		else		{
			rou0 = 1;
		}

		u = u + rou0;
	}



L_Bound[index] = l/numPoint;
U_Bound[index] = u/numPoint;

   }

}


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	//    [L,U] = estimation_gpu(x,y,z,sigma,bestbranch,bits_per_dim,dimension,num_thread,num_thread_per_block,num_block);

	double* x, *y, *z, *bestBranch;
	double bits_per_dim, dimension, num_thread;
	double sigma;

	x = mxGetPr(prhs[0]);
	y = mxGetPr(prhs[1]); 
	z = mxGetPr(prhs[2]); 
	sigma = *(mxGetPr(prhs[3]));
	bestBranch = mxGetPr(prhs[4]); 
	bits_per_dim = *(mxGetPr(prhs[5]));
	dimension = *(mxGetPr(prhs[6]));
	num_thread = *(mxGetPr(prhs[7]));

	int pN;	
	pN = mxGetN(prhs[0]); 

	double *L_Bound, *U_Bound;

	plhs[0] = mxCreateDoubleMatrix(1, (int)num_thread, mxREAL); //输出
	plhs[1] = mxCreateDoubleMatrix(1, (int)num_thread, mxREAL); //输出
	L_Bound = mxGetPr(plhs[0]);
	U_Bound = mxGetPr(plhs[1]);

	getLU(x, y, z, sigma, bestBranch, bits_per_dim, dimension, num_thread, pN, L_Bound, U_Bound);


//	for (int ii = 0;ii<=15;ii++){
//		printf("%f\t",L_Bound[ii]);
//		printf("%f\n",U_Bound[ii]);
//	}



}


double getmaxormin(double x, double paramax, double paramin, int flag){
    double result = 0;
    if (flag==1){     // 求ax的最大值

        if (x>=0)
        {result = paramax*x;}
        else
        {result = paramin*x;}
                
    }
    
    if (flag==-1){   //求ax的最小值
        if (x>=0)
        {result = paramin*x;}
        else
        {result = paramax*x;}
    }
        return result;

}





