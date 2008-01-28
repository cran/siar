// Usefile functions for the geese data set

#include"use.h"
#include<Rmath.h>
#include<R.h>

double GetLik(double *X,double mean,double sd,int datasize)
{
// This function gets the normal density for two objects with differing means and sds
    
double Ans=0.0;
int k;

for(k=0;k<datasize;k++)  {
    Ans += dnorm(X[k],mean,sd,1.0);
    //Rprintf("X[k] = %lf, mean=%lf, sd = %lf, Ans = %lf \n",X[k],mean,sd,Ans);
}

return(Ans);

}

//maximum function:

double maximum (double e, double f)
{
    if (e > f)
        return e;
    else
        return f;
}

// RW proposal functions
//rtruncn function:
double rtruncn (double a, double b)
{
    double A, B;
    double maxA, maxB, maxR, r2, r, th, u, x;
    int accept=0;
    
    A = atan(a);
    B = atan(b);
    
    maxA = exp(-pow(a,2)/4)/cos(A);
    maxB = exp(-pow(b,2)/4)/cos(B);
    maxR = maximum(maxA, maxB);
    if((a<1) && (b>-1))
        maxR = exp(-0.25)*sqrt(2.0);
    while (accept==0)
    {
        r2 = runif(0.0,1.0);
        r = sqrt(r2)*maxR;
        th = runif(A,B);
        u = r*cos(th);
        x = tan(th);
        if((pow(x,2)) < (log(u)*-4)) accept = 1;
    }
    return x;
}        


////////////////////////////////////////////////////////////////////////////////

//truncated normal function:
double truncatedwalk (double old, double sd, double low, double high)
{
    double lowlimold, upplimold, y, newvalue;
    lowlimold = (low - old)/sd;
    upplimold = (high - old)/sd;
    y = rtruncn(lowlimold, upplimold);
    newvalue = old + sd*y;
           
    return newvalue;
}

////////////////////////////////////////////////////////////////////////////////

//truncated normal ratio function:
double truncatedrat (double old, double sd, double low, double high, double newvalue)
{
    double lowlimold, upplimold, lowlimnew, upplimnew, plowold, puppold, plownew, puppnew, ratio;
    
    lowlimold = (low - old)/sd;
    upplimold = (high - old)/sd;
    lowlimnew = (low - newvalue)/sd;
    upplimnew = (high - newvalue)/sd;
    //plowold = normal_cdf(lowlimold);
    //puppold = normal_cdf(upplimold);
    //plownew = normal_cdf(lowlimnew);
    //puppnew = normal_cdf(upplimnew);
    plowold = pnorm(lowlimold,0,1,0,0);
    puppold = pnorm(upplimold,0,1,0,0);
    plownew = pnorm(lowlimnew,0,1,0,0);
    puppnew = pnorm(upplimnew,0,1,0,0);
    ratio = (puppold - plowold)/(puppnew - plownew);
    return ratio;        
}


double UpdateMCMC(double newloglik,double oldloglik,double newval,double oldval,double rat)
{
// Function to update MCMC when given log likelihoods

double u,mh;
u = runif(0.0,1.0);
mh = exp(newloglik - oldloglik)*rat;
if (u <= mh)
{
	return(newval);
} else 
{
	return(oldval);
}

}

float gammln(float xx)
// Returns the value log(Gamma(xx)) for xx>0
{
    double x,y,tmp,ser;
    static double cof[6] = {76.18009172947146,-86.50532032941677,
        24.01409824083091,-1.231739572450155,
        0.1208650973866179e-2,-0.5395239384953e-5};
    int j;

    y=x=xx;
    tmp=x+5.5;
    tmp -= (x+0.5)*log(tmp);
    ser=1.000000000190015;
    for(j=0;j<=5;j++) ser+=cof[j]/++y;
    return -tmp+log(2.5066282746310005*ser/x);
}

double logddirichlet(double *x,double *alpha,int len)
{

// This function calculates teh log dirichlet density
double logB=0.0,logdens,sumalpha=0.0,s=0.0;
int k;

for(k=0;k<len;k++) {
    s += (alpha[k]-1)*x[k];
    sumalpha += alpha[k];
    logB += gammln(alpha[k]);
}

logdens = s-logB-gammln(sumalpha);

return(logdens);

}
