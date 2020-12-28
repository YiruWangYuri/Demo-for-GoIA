#include "mex.h"
#include "math.h"
#include<stdio.h>
void getLU(double *x, double *y, double *z, double *newBranch, int M,double sigma,int pN, double *l, double *u)
{
    double maxaa(double x, double y);
    double minaa(double x, double y);
    int i =0,j=0;
    int Dim = 8;
    double rou,rou0;
    double centerBranch[4]={0,0,0,0};//?????

    double a_min,a_max,b_min,b_max,c_min,c_max,d_min,d_max,e_max,e_min,r_max,r_min,e0;
 
    for (i=0;i<M;i++)
    {
        l[i]=0;u[i]=0;rou=0;rou0=0;
            
        a_min=newBranch[0+i*Dim];
        a_max=newBranch[4+i*Dim];
        b_min=newBranch[1+i*Dim];
        b_max=newBranch[5+i*Dim];
        c_min=newBranch[2+i*Dim];
        c_max=newBranch[6+i*Dim];
        d_min=newBranch[3+i*Dim];
        d_max=newBranch[7+i*Dim];
                
//         printf("a_min %lf\n",a_min);
//         printf("a_max %lf\n",a_max);
//         printf("b_min %lf\n",b_min);
//         printf("b_max %lf\n",b_max);
//         printf("c_min %lf\n",c_min);
//         printf("c_max %lf\n",c_max);
//         printf("d_min %lf\n",d_min);
//         printf("d_max %lf\n",d_max);
        
        centerBranch[0]=0.5*(a_min+a_max);
        centerBranch[1]=0.5*(b_min+b_max);
        centerBranch[2]=0.5*(c_min+c_max);
        centerBranch[3]=0.5*(d_min+d_max);
        
//         printf("centerBranch0 %d\n",centerBranch[0]);
        
        for (j = 0;j<pN;j++)
        {
            e_max=a_max*x[j]+b_max*y[j]+c_max*z[j]+d_max;
            e_min=a_min*x[j]+b_min*y[j]+c_min*z[j]+d_min;
            if(e_max*e_min<0)
                {r_min = 0;
            
                e_max=fabs(e_max);
                e_min=fabs(e_min);
                r_max = maxaa(e_max,e_min);}
            else                
            {   e_max=fabs(e_max);
                e_min=fabs(e_min);
                r_min=minaa(e_max,e_min);
                r_max=maxaa(e_max,e_min);}
            
            if (r_min>sigma)
            {
                rou = 1;
            }
            else 
            {   
                rou = 0;
            }
             
        
            
            l[i] = l[i]+rou;
            
            e0=fabs(centerBranch[0]*x[j]+centerBranch[1]*y[j]+centerBranch[2]*z[j]+centerBranch[3]);
            if(e0<=sigma)
                rou0=0;
            else
                rou0=1;
            
            u[i]=u[i]+rou0;

            

        }
//         printf("%d\n",pN);
//        printf("ui %lf\n",u[i]);
        l[i] = l[i]/pN;
        u[i] = u[i]/pN;
    }
    
    
    

//     l[0] = x[0];l[1] = l[1];l[2] = l[2];l[3] = x[3];
//     u[0] = y[0];u[1] = y[1];u[2] = y[2];u[3] = y[3];

}

double maxaa(double x, double y){
    if (x>=y)
       return x;
    else
        return y;
}

double minaa(double x, double y){
    if (x<=y)
       return x;
    else
       return y;
}

void mexFunction(int nlhs,mxArray *plhs[], int nrhs,const mxArray *prhs[])
{
double *a, *b;
double *x, *y, *z, *newBranch;
int M,pN;double sigma;
x = mxGetPr(prhs[0]);
y = mxGetPr(prhs[1]);
z = mxGetPr(prhs[2]);
newBranch  = mxGetPr(prhs[3]);
M = *(mxGetPr(prhs[4]));
pN = mxGetN(prhs[0]);
sigma = *(mxGetPr(prhs[5]));

plhs[0] = mxCreateDoubleMatrix(1, M, mxREAL);
plhs[1] = mxCreateDoubleMatrix(1, M, mxREAL);
a = mxGetPr(plhs[0]);
b = mxGetPr(plhs[1]);


getLU(x, y,z,newBranch,M,sigma,pN,a,b);
//         printf("%d\n",newBranch[0]);
//         printf("%d\n",newBranch[1]);
//         printf("%d\n",newBranch[2]);
//         printf("%d\n",newBranch[3]);
//         printf("%d\n",newBranch[4]);
//         printf("%d\n",newBranch[5]);
//         printf("%d\n",newBranch[6]);
//         printf("%d\n",newBranch[7]);
}
