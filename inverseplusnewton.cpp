#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>
#include <cstdlib>
#include "f.h"
#define a0 2.50662823884
#define a1 -18.61500062529
#define a2 41.39119773534
#define a3 -25.44106049637
#define c0 0.3374754822726147
#define c1 0.9761690190917186
#define c2 0.1607979714918209
#define c3 0.0276438810333863
#define c4 0.0038405729373609
#define b0 -8.47351093090
#define b1 23.08336743743
#define b2 -21.06224101826
#define b3 3.13082909833
#define c5 0.0003951896511919
#define c6 0.0000321767881768
#define c7 0.0000002888167364
#define c8 0.0000003960315187
#define c 0.918938533204672








using namespace std;
#define pi 3.1415926

clock_t t;

double N(const double& x) {
    //    // double x;
    //    // x=30*get_uniform()-15;
    //
    double j,z,h;
    double a,b,y;
    j=floor(min(abs(x) + 0.5, 14));
    //  cout<<j<<endl;
    
    z = j;
    h =abs(x)-z;
    // if(x==15)
    // cout<<h<<endl;
    a =get_vvalue(j);
    b=z*a-1.0;
    //   cout<<b<<endl;
    
    double q=1.0;
    double s=0.0;
    s=a+h*b;
    for(double i=2.0;i<=(24-j);)
    {
        a=(a +z*b)/i;
        b=(b+z*a)/(i+1);
        q=q*h*h;
        // cout<<q;
        s=s+q*(a+h*b);
        i=i+2.0;
    }
    
    y=(s*exp(-0.5*x*x-c));
    if(x>0)
        y=1-y;
    return y;
    //
    
    
    
    //    double z=x;
    //    if (z > 6.0) { return 1.0; }; // this guards against overflow
    //    if (z < -6.0) { return 0.0; };
    //    double d1 = 0.31938153;
    //    double d2 = -0.356563782;
    //    double d3 = 1.781477937;
    //    double d4 = -1.821255978;
    //    double d5 = 1.330274429;
    //    double p = 0.2316419;
    //    double e2 = 0.3989423;
    //    double a=fabs(z);
    //    double t = 1.0/(1.0+a*p);
    //    double b = e2*exp((-z)*(z/2.0));
    //    double n = ((((d5*t+d4)*t+d3)*t+d2)*t+d1)*t;
    //    n = 1.0-b*n;
    //    if ( z < 0.0 ) n = 1.0 - n;
    //    return n;
    
    //cout<<y<<endl;
    
    
    //
    
    
};


double get_gaussian()

{
    double u=get_uniform();
   // cout<<u<<endl;
    double y=u-0.5;
    double r,x;
    if(abs(y)<0.42)
    {
        r=y*y;
        x=y*(((a3*r+a2)*r+a1)*r+a0)/((((b3*r+b2)*r+b1)*r+b0)*r+1);
        
    }
    else
    {
        r=u;
        if (y>0)
            r=1-u;
        r=log(-log(r));
        x=c0 +r*(c1 +r*(c2 +r*(c3 +r*(c4+r*(c5 +r*(c6 +r*(c7 +r*c8)))))));
        if (y < 0)
            x=-x;
        
    }
    x=x+(u-N(x))*exp(-0.5*x*x+c);
     //cout<<x<<endl;
    return x;
}

int main (int argc, char* argv[])
{
    double Si,r,sigma,T,K,St,temp,average,xbar,ybar;
    r=0.1215;
    sigma=0.4;
    Si=8;
    T=2.0;
    K=8;
    
    
    double putprice,sum;
    // sum=0.0;
    
    temp=pow(sigma,2.0);
    
    double interval_upperbound, interval_lowerbound,quantile,estimated_standard_error,width_of_confidence_interval;
   // quantile=1.96;
    double n=10000000.0;
    //width_of_confidence_interval=1.0;
    
    sum=0.0;
    double k=2.0;
    for(int i=1;i<=n;i++)
    {
        
        St=Si*exp ( (r-(temp/2.0) )*T+sigma*sqrt(T)*get_gaussian());
//        cout<<"ccccc"<<sigma*sqrt(T)*get_gaussian()<<endl;
//        cout<<"dddd"<<(r-(temp/2.0) )*T<<endl;
//        cout<<"eeee"<<exp ( (r-(temp/2.0) )*T+sigma*sqrt(T)*get_gaussian())<<endl;

        
        putprice=exp(-r*T)*(max(0.0,(K-St)));
        sum+=putprice;
        if(i==1)
        {
            xbar=putprice;
            ybar=pow(putprice, 2.0);
            continue;
        }
        
        xbar=(1-(1/k))*xbar+(1/k)*putprice;
        ybar=(1-(1/k))*ybar+(1/k)*pow(putprice,2.0);
        k=k+1.0;
        
    }
    estimated_standard_error=sqrt((1/(n-1))*(ybar-pow(xbar,2.0)));
    average=sum/(n);
    
    // estimated_standard_error=sigma/sqrt(n);
   // interval_lowerbound=average-quantile*estimated_standard_error;
   // interval_upperbound=average+quantile*estimated_standard_error;
   // width_of_confidence_interval=interval_upperbound-interval_lowerbound;
    
    
    
    cout<<"number of trials :"<<n<<endl;
    cout<<"European option put price : "<<average<<endl;
    cout<<"absolute errors: "<<1.6547-average<<endl;
    cout<<"estimated standard error :"<<estimated_standard_error<<endl;
   // cout<<"95% confidence interval : ["<<interval_lowerbound<<" , "<<interval_upperbound<<"]"<<endl;
   // cout<<"width of confidence interval :"<<width_of_confidence_interval<<endl;
    t = clock() - t;
    cout<<"time consumed :"<<(float)(t)/float(CLOCKS_PER_SEC)<<" seconds"<<endl;
   // cout<<"efficiency :"<<pow(estimated_standard_error,2.0)*(float)(t)/float(CLOCKS_PER_SEC)<<endl;
    return 0;
    
}
