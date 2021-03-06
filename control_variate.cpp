#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>
#include <cstdlib>
#include "Header.h"

using namespace std;
#define pi 3.1415926

clock_t t;

//double get_uniform()
//{
//    return (((double) random())/(pow(2.0, 31.0)-1.0));
//}

double get_gaussian()

{
    double R,V,Z;
    //R=-2.0*log(get_uniform());
    R=-2.0*log(Uniform01());
    //V=2*pi*get_uniform();
    V=2*pi*Uniform01();
    Z=(sqrt(R)*cos(V));
    
    
    return Z ;
}

double max(double a, double b) {
    return (b < a )? a:b;
}

int main (int argc, char* argv[])
{
    
    
    double Si,r,q,sigma,T,K,St,temp,average,xbar,ybar;
    r=0.5/100;
    sigma=0.3;
    Si=2000;
    T=1.0/12.0;
    K=2200;
    q=2.0/100;
    
    double callprice,sum,Xt,m;
    // sum=0.0;
    m=0;
    temp=pow(sigma,2.0);
    
    double interval_upperbound, interval_lowerbound,quantile,estimated_standard_error,width_of_confidence_interval;
    //quantile=1.96;
    double n=10000000.0;
    //width_of_confidence_interval=1.0;
    callprice=0;
    sum=0.0;
    double k=2.0;

    
//    double Stbar,xtbar,sum1,sum2;
//    sum1=0.0;
//    sum2=0.0;
//    double xi[1000];
//    double yi[1000];
//    double b1,b2,b,sum3,sum4;
//    for(int i=1;i<=1000;i++)
//    {
//        Xt=(r-q-(temp/2.0))*T+sigma*sqrt(T)*get_gaussian();
//        //m=get_gaussian();
//        //cout<<m<<endl;
//        St=Si*exp(Xt);
//        if(St>K)
//            
//            callprice=exp(-r*T)*St;
//        else
//            callprice=0;
//        xi[i]=St;
//        yi[i]=callprice;
//        sum1=St+sum1;
//        sum2+=callprice;
//    }
//    
//    for(int i=1;i<=1000;i++)
//    {
//        b1+=(xi[i]-sum1/1000)*(yi[i]-sum2/1000);
//        b2+=(xi[i]-sum1/1000)*(xi[i]-sum1/1000);
//    }
//    b=b1/b2;
//    cout<<b<<endl;

    double b=3.10622;
    for(int i=1;i<=n;i++)
    {
        Xt=(r-q-(temp/2.0))*T+sigma*sqrt(T)*get_gaussian();
        //m=get_gaussian();
        //cout<<m<<endl;
        St=Si*exp(Xt);
        if(St>K)
            
            callprice=exp(-r*T)*St+b*(Si*exp((r-q)*T)-St);
        else
            callprice=0+b*(Si*exp((r-q)*T)-St);
        
        //cout<<Xt<<endl;
        
        //cout<<St<<endl;
        //callprice=exp(-r*T)*(max((K-St),0.0));
        sum+=callprice;
        //sum+=St;
        
        if(i==1)
        {
            xbar=callprice;
            ybar=pow(callprice, 2.0);
            continue;
        }
        
        xbar=(1-(1/k))*xbar+(1/k)*callprice;
        ybar=(1-(1/k))*ybar+(1/k)*pow(callprice,2.0);
        k=k+1.0;
        
    }
    estimated_standard_error=sqrt((1/(n-1))*(ybar-pow(xbar,2.0)));
    average=sum/(n);
    
    // estimated_standard_error=sigma/sqrt(n);
    // interval_lowerbound=average-quantile*estimated_standard_error;
    //interval_upperbound=average+quantile*estimated_standard_error;
    // width_of_confidence_interval=interval_upperbound-interval_lowerbound;
    
    cout<<"number of trials :"<<n<<endl;
    cout<<"asset or nothing option call price : "<<average<<endl;
    //cout<<"absolute errors: "<<1463.06-average<<endl;
    cout<<"absolute errors: "<<283.391-average<<endl;

    cout<<"estimated standard error :"<<estimated_standard_error<<endl;
    
    // cout<<"95% confidence interval : ["<<interval_lowerbound<<" , "<<interval_upperbound<<"]"<<endl;
    // cout<<"width of confidence interval :"<<width_of_confidence_interval<<endl;
    t = clock() - t;
    cout<<"time consumed :"<<(float)(t)/float(CLOCKS_PER_SEC)<<" seconds"<<endl;
    cout<<"efficiency :"<<pow(estimated_standard_error,2.0)*(float)(t)/float(CLOCKS_PER_SEC)<<endl;
    cout<<m;
    return 0;
    
}
